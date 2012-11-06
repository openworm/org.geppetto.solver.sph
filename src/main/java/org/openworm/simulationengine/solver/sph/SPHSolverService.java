package org.openworm.simulationengine.solver.sph;

import static java.lang.System.out;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.bridj.Pointer;
import org.openworm.simulationengine.core.constants.PhysicsConstants;
import org.openworm.simulationengine.core.model.IModel;
import org.openworm.simulationengine.core.simulation.ITimeConfiguration;
import org.openworm.simulationengine.core.solver.ISolver;
import org.openworm.simulationengine.model.sph.SPHParticle;
import org.openworm.simulationengine.model.sph.common.SPHConstants;
import org.openworm.simulationengine.model.sph.x.SPHModelX;
import org.openworm.simulationengine.model.sph.x.SPHParticleX;
import org.openworm.simulationengine.model.sph.x.Vector3DX;
import org.springframework.stereotype.Service;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.opencl.CLPlatform.DeviceFeature;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.util.IOUtils;

@Service
public class SPHSolverService implements ISolver {
	
	
	private static Log logger = LogFactory.getLog(SPHSolverService.class);
	
	private CLContext _context;
	public CLQueue _queue;
	private CLProgram _program;
	private CLDevice _device;
	public CLBuffer<Float> _acceleration;
	public CLBuffer<Integer> _gridCellIndex;
	public CLBuffer<Integer> _gridCellIndexFixedUp;
	public CLBuffer<Float> _neighborMap;
	public CLBuffer<Integer> _particleIndex;
	public CLBuffer<Float> _position;
	public CLBuffer<Float> _pressure;
	public CLBuffer<Float> _rho;
	public CLBuffer<Float> _rhoInv;
	public CLBuffer<Float> _sortedPosition;
	public CLBuffer<Float> _sortedVelocity;
	public CLBuffer<Float> _velocity;
	
	public Pointer<Float> _accelerationPtr;
	public Pointer<Integer> _gridCellIndexPtr;
	public Pointer<Integer> _gridCellIndexFixedUpPtr;
	public Pointer<Float> _neighborMapPtr;
	public Pointer<Integer> _particleIndexPtr;
	public Pointer<Float> _positionPtr;
	public Pointer<Float> _positionPtrbuff;
	public Pointer<Float> _pressurePtr;
	public Pointer<Float> _rhoPtr;
	public Pointer<Float> _rhoInvPtr;
	public Pointer<Float> _sortedPositionPtr;
	public Pointer<Float> _sortedVelocityPtr;
	public Pointer<Float> _velocityPtr;
	
	private CLKernel _clearBuffers;
	private CLKernel _computeAcceleration;
	private CLKernel _computeDensityPressure;
	private CLKernel _findNeighbors;
	private CLKernel _hashParticles;
	private CLKernel _indexPostPass;
	private CLKernel _indexx;
	private CLKernel _integrate;
	private CLKernel _sortPostPass;
	
	public int _gridCellsX;
	public int _gridCellsY;
	public int _gridCellsZ;
	public int _gridCellCount;
	public int _particleCount;
	
	public static Random RandomGenerator = new Random();
	
	public SPHSolverService() throws Exception{
		this.onceOffInit();
	}
		
	private void onceOffInit() throws IOException  
	{
		_context = JavaCL.createBestContext(DeviceFeature.CPU);
		
		out.println("created "+ _context);
		// an array with available devices
		CLDevice[] devices = _context.getDevices();

		for(int i=0; i<devices.length; i++)
		{
			out.println("device - " + i + ": " + devices[i]);
		}	

		// have a look at the output and select a device
		_device = devices[0];
		out.println("Version " + _device.getOpenCLVersion());
		out.println("Version " + _device.getDriverVersion());
		out.println("using "+ _device);
		out.println("max workgroup size: " + _device.getMaxWorkGroupSize());
		out.println("max workitems size: " + _device.getMaxWorkItemSizes()[0]);
		
		// create command queue on selected device.
		_queue = _context.createDefaultQueue();//device.createCommandQueue();
		
		// load sources, create and build program
		String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/resource/sphFluidDemo.cl"));
		_program = _context.createProgram(src);
		
		/*kernels*/
		_clearBuffers = _program.createKernel("clearBuffers");
		_computeAcceleration = _program.createKernel("computeAcceleration");
		_computeDensityPressure = _program.createKernel("computeDensityPressure");
		_findNeighbors = _program.createKernel("findNeighbors");
		_hashParticles = _program.createKernel("hashParticles");
		_indexPostPass = _program.createKernel("indexPostPass");
		_indexx = _program.createKernel("indexx");
		_integrate = _program.createKernel("integrate");
		_sortPostPass = _program.createKernel("sortPostPass");
	}
	
	private void allocateBuffers(){
		// input buffers declarations
		_accelerationPtr = Pointer.allocateFloats(_particleCount * 4);
		_gridCellIndexPtr = Pointer.allocateInts((_gridCellCount + 1));
		_gridCellIndexFixedUpPtr = Pointer.allocateInts((_gridCellCount + 1));
		_neighborMapPtr = Pointer.allocateFloats(_particleCount * SPHConstants.NEIGHBOR_COUNT * 2);
		_particleIndexPtr = Pointer.allocateInts(_particleCount * 2);
		_positionPtr = Pointer.allocateFloats(_particleCount * 4);
		_pressurePtr = Pointer.allocateFloats(_particleCount * 1);
		_rhoPtr = Pointer.allocateFloats(_particleCount * 1);
		_rhoInvPtr = Pointer.allocateFloats(_particleCount * 1);
		_sortedPositionPtr = Pointer.allocateFloats(_particleCount * 4);
		_sortedVelocityPtr = Pointer.allocateFloats(_particleCount * 4);
		_velocityPtr = Pointer.allocateFloats(_particleCount * 4);
		_positionPtrbuff = Pointer.allocateFloats(_particleCount * 4);
		
		// alternative buffer defining
		_acceleration = _context.createBuffer(Usage.InputOutput,_accelerationPtr,false);
		_gridCellIndex = _context.createBuffer(Usage.InputOutput,_gridCellIndexPtr, false);
		_gridCellIndexFixedUp = _context.createIntBuffer(Usage.InputOutput,_gridCellIndexFixedUpPtr, false);
		_neighborMap = _context.createBuffer(Usage.InputOutput,_neighborMapPtr, false);
		_particleIndex = _context.createBuffer(Usage.InputOutput,_particleIndexPtr ,false);
		_position = _context.createBuffer(Usage.InputOutput, _positionPtr, false);
		_pressure = _context.createBuffer(Usage.InputOutput,_pressurePtr,false);
		_rho = _context.createBuffer(Usage.InputOutput,_rhoPtr,false);
		_rhoInv = _context.createBuffer(Usage.InputOutput,_rhoInvPtr,false);
		_sortedPosition = _context.createBuffer(Usage.InputOutput,_sortedPositionPtr,false);
		_sortedVelocity = _context.createBuffer(Usage.InputOutput, _sortedVelocityPtr,false);
		_velocity = _context.createBuffer(Usage.InputOutput,_velocityPtr,false);
		_queue.finish();
	}
	
	public void setModels(List<IModel> models){
		// TODO: generalize this for an arbitrary number of models instead of just one
		if(!(models == null || models.size() ==0))
		{
			SPHModelX mod = (SPHModelX) models.get(0);
			
			_particleCount = mod.getNumberOfParticals();
			
			// set grid dimensions
			_gridCellsX = mod.getCellX();
			_gridCellsY = mod.getCellY();
			_gridCellsZ = mod.getCellZ();
			_gridCellCount = _gridCellsX * _gridCellsY * _gridCellsZ;
			
			// allocate buffers - requires global dimensions of the grid
			this.allocateBuffers();
			
			int index = 0;
			
			for(int i = 0;i< _particleCount;i++){
				if(i != 0)
				{
					index = index + 4;
				}
				
				Vector3DX positionVector = (Vector3DX) mod.getParticles().get(i).getPositionVector();
				Vector3DX velocityVector = (Vector3DX) mod.getParticles().get(i).getVelocityVector();
				
				// buffer population
				_positionPtr.set(index,positionVector.getX());
				_positionPtr.set(index + 1,positionVector.getY());
				_positionPtr.set(index + 2,positionVector.getZ());
				_positionPtr.set(index + 3,positionVector.getP());
				_velocityPtr.set(index,velocityVector.getX());
				_velocityPtr.set(index + 1,velocityVector.getY());
				_velocityPtr.set(index + 2,velocityVector.getZ());
				_velocityPtr.set(index + 3,velocityVector.getP());
			}
		}
		else
		{
			throw new IllegalArgumentException("SPHSolverService:setModels - invalid models");
		}
	}
	
	public List<IModel> getModels(){
		List<IModel> models = new ArrayList<IModel>();
		
		SPHModelX mod = new SPHModelX(_gridCellsX, _gridCellsY, _gridCellsZ);
		
		int index = 0;
		for(int i = 0;i< _particleCount;i++){
			if(i != 0)
			{
				index = index + 4;
			}
			
			Vector3DX positionVector = new Vector3DX();
			Vector3DX velocityVector = new Vector3DX();
			
			// buffer population
			positionVector.setX(_positionPtr.get(index));
			positionVector.setY(_positionPtr.get(index + 1));
			positionVector.setZ(_positionPtr.get(index + 2));
			positionVector.setP(_positionPtr.get(index + 3));
			
			velocityVector.setX(_velocityPtr.get(index));
			velocityVector.setY(_velocityPtr.get(index + 1));
			velocityVector.setZ(_velocityPtr.get(index + 2));
			velocityVector.setP(_velocityPtr.get(index + 3));
			
			SPHParticle particle = new SPHParticleX(positionVector, velocityVector, 1);
			((SPHParticleX)particle).setId(mod.getId()+i);
			mod.getParticles().add(particle);
		}
		
		models.add(mod);
		
		return models;
	}
	
	public List<List<IModel>> solve(List<IModel> models, ITimeConfiguration timeConfiguration)
	{
		// TODO: extend this to use time configuration to do multiple steps in one go
		logger.info("SPH solver start");
		// 1. populate buffers from list of models
		setModels(models);
		
		// 2. call this.step();
		step();
		
		// 3. retrieve values from buffers and populate returned models
		List<List<IModel>> modelsList = new ArrayList<List<IModel>> ();
		modelsList.add(this.getModels());
		logger.info("SPH solver end");
		return modelsList;
	}
	
	public void cleanContext(){
		this._context.release();
	}
	
	public int runClearBuffers(){
		_clearBuffers.setArg(0, _neighborMap);
		_clearBuffers.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runComputeAcceleration(){
		_computeAcceleration.setArg( 0, _neighborMap );
		_computeAcceleration.setArg( 1, _pressure );
		_computeAcceleration.setArg( 2, _rho );
		_computeAcceleration.setArg( 3, _rhoInv );
		_computeAcceleration.setArg( 4, _sortedPosition );
		_computeAcceleration.setArg( 5, _sortedVelocity );
		_computeAcceleration.setArg( 6, PhysicsConstants.CFLLimit );
		_computeAcceleration.setArg( 7, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT );
		_computeAcceleration.setArg( 8, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		_computeAcceleration.setArg( 9, PhysicsConstants.H );
		_computeAcceleration.setArg( 10, PhysicsConstants.MASS );
		_computeAcceleration.setArg( 11, PhysicsConstants.MU );
		_computeAcceleration.setArg( 12, PhysicsConstants.SIMULATION_SCALE );
		_computeAcceleration.setArg( 13, _acceleration );
		_computeAcceleration.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runComputeDensityPressure(){
		_computeDensityPressure.setArg( 0, _neighborMap );
		_computeDensityPressure.setArg( 1, PhysicsConstants.W_POLY_6_COEFFICIENT );
		_computeDensityPressure.setArg( 2, PhysicsConstants.H );
		_computeDensityPressure.setArg( 3, PhysicsConstants.MASS );
		_computeDensityPressure.setArg( 4, PhysicsConstants.RHO0 );
		_computeDensityPressure.setArg( 5, PhysicsConstants.SIMULATION_SCALE );
		_computeDensityPressure.setArg( 6, PhysicsConstants.STIFFNESS );
		_computeDensityPressure.setArg( 7, _pressure );
		_computeDensityPressure.setArg( 8, _rho );
		_computeDensityPressure.setArg( 9, _rhoInv );
		_computeDensityPressure.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runFindNeighbors(){
		_findNeighbors.setArg( 0, _gridCellIndexFixedUp );
		_findNeighbors.setArg( 1, _sortedPosition );
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_findNeighbors.setArg( 2, _gridCellCount );
		_findNeighbors.setArg( 3, _gridCellsX );
		_findNeighbors.setArg( 4, _gridCellsY );
		_findNeighbors.setArg( 5, _gridCellsZ );
		_findNeighbors.setArg( 6, PhysicsConstants.H );
		_findNeighbors.setArg( 7, PhysicsConstants.HASH_GRID_CELL_SIZE );
		_findNeighbors.setArg( 8, PhysicsConstants.HASH_GRID_CELL_SIZE_INV );
		_findNeighbors.setArg( 9, PhysicsConstants.SIMULATION_SCALE );
		_findNeighbors.setArg( 10, SPHConstants.XMIN_F );
		_findNeighbors.setArg( 11, SPHConstants.YMIN_F );
		_findNeighbors.setArg( 12, SPHConstants.ZMIN_F );
		_findNeighbors.setArg( 13, _neighborMap );
		_findNeighbors.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runHashParticles(){
		// Stage HashParticles
		_hashParticles.setArg( 0, _position );
		_hashParticles.setArg( 1, _gridCellsX );
		_hashParticles.setArg( 2, _gridCellsY );
		_hashParticles.setArg( 3, _gridCellsZ );
		_hashParticles.setArg( 4, PhysicsConstants.HASH_GRID_CELL_SIZE_INV );
		_hashParticles.setArg( 5, SPHConstants.XMIN_F );
		_hashParticles.setArg( 6, SPHConstants.YMIN_F );
		_hashParticles.setArg( 7, SPHConstants.ZMIN_F );
		_hashParticles.setArg( 8, _particleIndex );
		_hashParticles.setArg( 9, _particleCount );
		_hashParticles.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runIndexPostPass(){
		_indexPostPass.setArg( 0, _gridCellIndex );
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_indexPostPass.setArg( 1, _gridCellCount );
		_indexPostPass.setArg( 2, _gridCellIndexFixedUp );
		int gridCellCountRoundedUp = ((( _gridCellCount - 1 ) / 256 ) + 1 ) * 256;
		_indexPostPass.enqueueNDRange(_queue, new int[] {gridCellCountRoundedUp});
		return 0;
	}
	
	public int runIndexx(){
		// Stage Indexx
		_indexx.setArg( 0, _particleIndex );
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_indexx.setArg( 1, _gridCellCount );
		_indexx.setArg( 2, _gridCellIndex );
		_indexx.setArg( 3, _particleCount );
		int gridCellCountRoundedUp = ((( _gridCellCount - 1 ) / 256 ) + 1 ) * 256;
		_indexx.enqueueNDRange(_queue, new int[] {gridCellCountRoundedUp});
		return 0;
	}
	
	public int runIntegrate(){
		// Stage Integrate
		_integrate.setArg( 0, _acceleration );
		_integrate.setArg( 1, _sortedPosition );
		_integrate.setArg( 2, _sortedVelocity );
		_integrate.setArg( 3, PhysicsConstants.GRAVITY_X );
		_integrate.setArg( 4, PhysicsConstants.GRAVITY_Y );
		_integrate.setArg( 5, PhysicsConstants.GRAVITY_Z );
		_integrate.setArg( 6, PhysicsConstants.SIMULATION_SCALE_INV );
		_integrate.setArg( 7, PhysicsConstants.TIME_STEP );
		_integrate.setArg( 8, SPHConstants.XMIN_F );
		_integrate.setArg( 9, SPHConstants.XMAX_F );
		_integrate.setArg( 10, SPHConstants.YMIN_F );
		_integrate.setArg( 11, SPHConstants.YMAX_F );
		_integrate.setArg( 12, SPHConstants.ZMIN_F );
		_integrate.setArg( 13, SPHConstants.ZMAX_F );
		_integrate.setArg( 14, PhysicsConstants.DAMPING );
		_integrate.setArg( 15, _position );
		_integrate.setArg( 16, _velocity );
		_integrate.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runSortPostPass(){
		// Stage SortPostPass
		_sortPostPass.setArg( 0, _particleIndex );
		_sortPostPass.setArg( 1, _position );
		_sortPostPass.setArg( 2, _velocity );
		_sortPostPass.setArg( 3, _sortedPosition );
		_sortPostPass.setArg( 4, _sortedVelocity );
		_sortPostPass.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int runSort(){
		//this version work with qsort
		int index = 0;
		List<int[]> particleIndex = new ArrayList<int[]>();
		Pointer<Integer> particleInd = _particleIndex.read(_queue);
		_queue.finish();
		for(int i = 0; i < _particleCount * 2;i+=2){
			int[] element = {particleInd.get(i), particleInd.get(i+1)};
			particleIndex.add(element);
		}
		Collections.sort(particleIndex, new MyCompare());
		for(int i = 0; i< particleIndex.size();i++){
			for(int j=0;j<2;j++){
				particleInd.set(index,particleIndex.get(i)[j]);
				index++;
			}
		}
		_particleIndex.write(_queue, particleInd, false);
		_queue.finish();
		return 0;
	}
	
	class MyCompare implements Comparator<int[]>{
		public int compare(int[] o1, int[] o2) {
			if( o1[0] < o2[0] ) return -1;
			if( o1[0] > o2[0]) return +1;
			return 0;
		}
	}
	
	private void step(){
		logger.info("SPH run clear buffer");
		runClearBuffers();
		logger.info("SPH run hash particles");
		runHashParticles();
		logger.info("SPH run sort");
		runSort();
		logger.info("SPH run sort post pass");
		runSortPostPass();
		logger.info("SPH run index");
		runIndexx();
		logger.info("SPH run index post pass");
		runIndexPostPass();
		logger.info("SPH run find neighbors");
		runFindNeighbors();
		logger.info("SPH runComputeDensityPressur");
		runComputeDensityPressure();
		logger.info("SPH runComputeAcceleration");
		runComputeAcceleration();
		logger.info("SPH runIntegrate");
		runIntegrate();
		logger.info("SPH position read");
		logger.info("Position element size: " + _position.getElementSize());
		logger.info("Position element count: " + _position.getElementCount());		
		logger.info("Queue: "+ _queue);
		_positionPtr = _position.read(_queue);
		logger.info("SPH finish queue");
		_queue.finish();
		logger.info("SPH step done");
	}
	
	public void finishQueue() {
		_queue.finish();
	}
}
;