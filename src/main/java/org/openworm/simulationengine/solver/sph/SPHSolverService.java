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
	public CLBuffer<Integer> _particleIndexBack;
	public CLBuffer<Float> _position;
	public CLBuffer<Float> _pressure;
	public CLBuffer<Float> _rho;
	public CLBuffer<Float> _sortedPosition;
	public CLBuffer<Float> _sortedVelocity;
	public CLBuffer<Float> _velocity;
	public CLBuffer<Float> _elasticConnectionsData;
	
	public Pointer<Float> _accelerationPtr;
	public Pointer<Integer> _gridCellIndexPtr;
	public Pointer<Integer> _gridCellIndexFixedUpPtr;
	public Pointer<Float> _neighborMapPtr;
	public Pointer<Integer> _particleIndexPtr;
	public Pointer<Integer> _particleIndexBackPtr;
	public Pointer<Float> _positionPtr;
	public Pointer<Float> _pressurePtr;
	public Pointer<Float> _rhoPtr;
	public Pointer<Float> _sortedPositionPtr;
	public Pointer<Float> _sortedVelocityPtr;
	public Pointer<Float> _velocityPtr;
	public Pointer<Float> _elasticConnectionsDataPtr;
	
	/*
	 * Kernel declarations
	 */
	private CLKernel _clearBuffers;
	private CLKernel _findNeighbors;
	private CLKernel _hashParticles;
	private CLKernel _indexPostPass;
	private CLKernel _indexx;
	private CLKernel _sortPostPass;
	
	// additional kernels for PCISPH
	private CLKernel _pcisph_computeDensity;
	private CLKernel _pcisph_computeForcesAndInitPressure;
	private CLKernel _pcisph_integrate;
	private CLKernel _pcisph_predictPositions;
	private CLKernel _pcisph_predictDensity;
	private CLKernel _pcisph_correctPressure;
	private CLKernel _pcisph_computePressureForceAcceleration;
	private CLKernel _pcisph_computeElasticForces;
	
	public int _gridCellsX;
	public int _gridCellsY;
	public int _gridCellsZ;
	public int _gridCellCount;
	public int _particleCount;
	public int _numOfLiquidP;
	public int _numOfElasticP;
	public int _numOfBoundaryP;
	
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
		String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/resource/sphFluid.cl"));
		_program = _context.createProgram(src);
		
		// kernels
		_clearBuffers = _program.createKernel("clearBuffers");
		_findNeighbors = _program.createKernel("findNeighbors");
		_hashParticles = _program.createKernel("hashParticles");
		// PORTING-NOTE: indexPostPass is gone from kernels in the latest version,logic moved from opencl to host code (this class)
		_indexPostPass = _program.createKernel("indexPostPass");
		_indexx = _program.createKernel("indexx");
		_sortPostPass = _program.createKernel("sortPostPass");
		
		// PCI-SPH specific 
		_pcisph_computeForcesAndInitPressure  = _program.createKernel("pcisph_computeForcesAndInitPressure");
		_pcisph_integrate  = _program.createKernel("pcisph_integrate");
		_pcisph_predictPositions  = _program.createKernel("pcisph_predictPositions");
		_pcisph_predictDensity  = _program.createKernel("pcisph_predictDensity");
		_pcisph_correctPressure  = _program.createKernel("pcisph_correctPressure");
		_pcisph_computePressureForceAcceleration  = _program.createKernel("pcisph_computePressureForceAcceleration");
		_pcisph_computeDensity  = _program.createKernel("pcisph_computeDensity");
		_pcisph_computeElasticForces  = _program.createKernel("pcisph_computeElasticForces");
	}
	
	private void allocateBuffers(){
		// input buffers declarations
		_accelerationPtr = Pointer.allocateFloats(_particleCount * 4 * 2);
		_gridCellIndexPtr = Pointer.allocateInts((_gridCellCount + 1));
		_gridCellIndexFixedUpPtr = Pointer.allocateInts((_gridCellCount + 1));
		_neighborMapPtr = Pointer.allocateFloats(_particleCount * SPHConstants.NEIGHBOR_COUNT * 2);
		_particleIndexPtr = Pointer.allocateInts(_particleCount * 2);
		_particleIndexBackPtr = Pointer.allocateInts(_particleCount);
		_positionPtr = Pointer.allocateFloats(_particleCount * 4);
		_pressurePtr = Pointer.allocateFloats(_particleCount);
		_rhoPtr = Pointer.allocateFloats(_particleCount * 2);
		_sortedPositionPtr = Pointer.allocateFloats(_particleCount * 4 * 2);
		_sortedVelocityPtr = Pointer.allocateFloats(_particleCount * 4 * 2);
		_velocityPtr = Pointer.allocateFloats(_particleCount * 4);
		
		// alternative buffer defining
		_acceleration = _context.createBuffer(Usage.InputOutput,_accelerationPtr,false);
		_gridCellIndex = _context.createBuffer(Usage.InputOutput,_gridCellIndexPtr, false);
		_gridCellIndexFixedUp = _context.createIntBuffer(Usage.InputOutput,_gridCellIndexFixedUpPtr, false);
		_neighborMap = _context.createBuffer(Usage.InputOutput,_neighborMapPtr, false);
		_particleIndex = _context.createBuffer(Usage.InputOutput,_particleIndexPtr ,false);
		_particleIndexBack = _context.createBuffer(Usage.InputOutput,_particleIndexBackPtr ,false);
		_position = _context.createBuffer(Usage.InputOutput, _positionPtr, false);
		_pressure = _context.createBuffer(Usage.InputOutput,_pressurePtr,false);
		_rho = _context.createBuffer(Usage.InputOutput,_rhoPtr,false);
		_sortedPosition = _context.createBuffer(Usage.InputOutput,_sortedPositionPtr,false);
		_sortedVelocity = _context.createBuffer(Usage.InputOutput, _sortedVelocityPtr,false);
		_velocity = _context.createBuffer(Usage.InputOutput,_velocityPtr,false);
		
		// init elastic connections buffer if we have any
		if(_numOfElasticP > 0){
			_elasticConnectionsDataPtr = Pointer.allocateFloats(_numOfElasticP * SPHConstants.NEIGHBOR_COUNT * 4);
			_elasticConnectionsData = _context.createBuffer(Usage.InputOutput,_elasticConnectionsDataPtr,false);
		}
		
		_queue.finish();
	}
	
	public void setModels(List<IModel> models) {
		// TODO: generalize this for an arbitrary number of models instead of just one
		if(!(models == null || models.size() ==0))
		{
			SPHModelX mod = (SPHModelX) models.get(0);
			
			_particleCount = mod.getNumberOfParticals();
			// PORTING-TODO: populate elastic connection buffers
			_numOfElasticP = 0;
			_numOfLiquidP = 0;
			_numOfBoundaryP = 0;
			
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
				_positionPtr.set(index, round(positionVector.getX(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_positionPtr.set(index + 1, round(positionVector.getY(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_positionPtr.set(index + 2, round(positionVector.getZ(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_positionPtr.set(index + 3, round(positionVector.getP(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_velocityPtr.set(index, round(velocityVector.getX(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_velocityPtr.set(index + 1, round(velocityVector.getY(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_velocityPtr.set(index + 2, round(velocityVector.getZ(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				_velocityPtr.set(index + 3, round(velocityVector.getP(), SPHConstants.DECIMAL_ROUNDING_FACTOR));
				
				// particle counts
				if (positionVector.getP() == SPHConstants.BOUNDARY_TYPE) {
					_numOfBoundaryP++;
				}
				else if (positionVector.getP() == SPHConstants.ELASTIC_TYPE) {
					_numOfElasticP++;
				}
				else if (positionVector.getP() == SPHConstants.LIQUID_TYPE) {
					_numOfLiquidP++;
				}
			}
			
			// check that counts are fine
			if(_particleCount != (_numOfBoundaryP + _numOfElasticP + _numOfLiquidP)){
				throw new IllegalArgumentException("SPHSolverService:setModels - particle counts do not add up");
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
		_findNeighbors.setArg( 10, SPHConstants.XMIN );
		_findNeighbors.setArg( 11, SPHConstants.YMIN );
		_findNeighbors.setArg( 12, SPHConstants.ZMIN );
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
		_hashParticles.setArg( 5, SPHConstants.XMIN );
		_hashParticles.setArg( 6, SPHConstants.YMIN );
		_hashParticles.setArg( 7, SPHConstants.ZMIN );
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
	
	public int runSortPostPass(){
		// Stage SortPostPass
		_sortPostPass.setArg( 0, _particleIndex );
		_sortPostPass.setArg(1, _particleIndexBack );
		_sortPostPass.setArg( 2, _position );
		_sortPostPass.setArg( 3, _velocity );
		_sortPostPass.setArg( 4, _sortedPosition );
		_sortPostPass.setArg( 5, _sortedVelocity );
		_sortPostPass.setArg( 6, _particleCount );
		_sortPostPass.enqueueNDRange(_queue, new int[] {_particleCount});
		return 0;
	}
	
	public int run_pcisph_computeDensity(){
		// Stage ComputeDensityPressure
		_pcisph_computeDensity.setArg( 0, _neighborMap );
		_pcisph_computeDensity.setArg( 1, PhysicsConstants.W_POLY_6_COEFFICIENT );
		_pcisph_computeDensity.setArg( 2, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		_pcisph_computeDensity.setArg( 3, PhysicsConstants.H );
		_pcisph_computeDensity.setArg( 4, PhysicsConstants.MASS );
		_pcisph_computeDensity.setArg( 5, PhysicsConstants.RHO0 );
		_pcisph_computeDensity.setArg( 6, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_computeDensity.setArg( 7, PhysicsConstants.STIFFNESS );
		_pcisph_computeDensity.setArg( 8, _sortedPosition );
		_pcisph_computeDensity.setArg( 9, _pressure );
		_pcisph_computeDensity.setArg(10, _rho );
		_pcisph_computeDensity.setArg(11, _particleIndexBack );
		_pcisph_computeDensity.setArg(12, PhysicsConstants.DELTA ); // calculated from constants
		_pcisph_computeDensity.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_computeForcesAndInitPressure(){
		_pcisph_computeForcesAndInitPressure.setArg( 0, _neighborMap );
		_pcisph_computeForcesAndInitPressure.setArg( 1, _rho );
		_pcisph_computeForcesAndInitPressure.setArg( 2, _pressure );
		_pcisph_computeForcesAndInitPressure.setArg( 3, _sortedPosition );
		_pcisph_computeForcesAndInitPressure.setArg( 4, _sortedVelocity );
		_pcisph_computeForcesAndInitPressure.setArg( 5, _acceleration );
		_pcisph_computeForcesAndInitPressure.setArg( 6, _particleIndexBack );
		_pcisph_computeForcesAndInitPressure.setArg( 7, PhysicsConstants.W_POLY_6_COEFFICIENT );
		_pcisph_computeForcesAndInitPressure.setArg( 8, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT );
		_pcisph_computeForcesAndInitPressure.setArg( 9, PhysicsConstants.H );
		_pcisph_computeForcesAndInitPressure.setArg(10, PhysicsConstants.MASS );
		_pcisph_computeForcesAndInitPressure.setArg(11, PhysicsConstants.MU );
		_pcisph_computeForcesAndInitPressure.setArg(12, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_computeForcesAndInitPressure.setArg(13, PhysicsConstants.GRAVITY_X );
		_pcisph_computeForcesAndInitPressure.setArg(14, PhysicsConstants.GRAVITY_Y );
		_pcisph_computeForcesAndInitPressure.setArg(15, PhysicsConstants.GRAVITY_Z );
		_pcisph_computeForcesAndInitPressure.setArg(16, _position );
		_pcisph_computeForcesAndInitPressure.setArg(17, _particleIndex );
		_pcisph_computeForcesAndInitPressure.setArg(18, _particleCount );
		_pcisph_computeForcesAndInitPressure.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_computeElasticForces(){
		_pcisph_computeElasticForces.setArg( 0, _neighborMap );
		_pcisph_computeElasticForces.setArg( 1, _sortedPosition );
		_pcisph_computeElasticForces.setArg( 2, _sortedVelocity );
		_pcisph_computeElasticForces.setArg( 3, _acceleration );
		_pcisph_computeElasticForces.setArg( 4, _particleIndexBack );
		_pcisph_computeElasticForces.setArg( 5, PhysicsConstants.H );
		_pcisph_computeElasticForces.setArg( 6, PhysicsConstants.MASS );
		_pcisph_computeElasticForces.setArg( 7, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_computeElasticForces.setArg( 8, _numOfElasticP);
		_pcisph_computeElasticForces.setArg( 9, _elasticConnectionsData );
		_pcisph_computeElasticForces.setArg( 10, _numOfBoundaryP );

		int numOfElasticPRoundedUp = ((( _numOfElasticP - 1 ) / 256 ) + 1 ) * 256;

		_pcisph_computeElasticForces.enqueueNDRange(_queue, new int[] {numOfElasticPRoundedUp});
		
		return 0;
	}
	
	public int run_pcisph_predictPositions(){
		_pcisph_predictPositions.setArg( 0, _acceleration );
		_pcisph_predictPositions.setArg( 1, _sortedPosition );
		_pcisph_predictPositions.setArg( 2, _sortedVelocity );
		_pcisph_predictPositions.setArg( 3, _particleIndex );
		_pcisph_predictPositions.setArg( 4, _particleIndexBack );
		_pcisph_predictPositions.setArg( 5, PhysicsConstants.GRAVITY_X );
		_pcisph_predictPositions.setArg( 6, PhysicsConstants.GRAVITY_Y );
		_pcisph_predictPositions.setArg( 7, PhysicsConstants.GRAVITY_Z );
		_pcisph_predictPositions.setArg( 8, PhysicsConstants.SIMULATION_SCALE_INV );
		_pcisph_predictPositions.setArg( 9, PhysicsConstants.TIME_STEP );
		_pcisph_predictPositions.setArg( 10, SPHConstants.XMIN );
		_pcisph_predictPositions.setArg( 11, SPHConstants.XMAX );
		_pcisph_predictPositions.setArg( 12, SPHConstants.YMIN );
		_pcisph_predictPositions.setArg( 13, SPHConstants.YMAX );
		_pcisph_predictPositions.setArg( 14, SPHConstants.ZMIN );
		_pcisph_predictPositions.setArg( 15, SPHConstants.ZMAX );
		_pcisph_predictPositions.setArg( 16, PhysicsConstants.DAMPING );
		_pcisph_predictPositions.setArg( 17, _position );
		_pcisph_predictPositions.setArg( 18, _velocity );
		_pcisph_predictPositions.setArg( 19, PhysicsConstants.R0 );
		_pcisph_predictPositions.setArg( 20, _neighborMap );
		_pcisph_predictPositions.setArg( 21, _particleCount );
		_pcisph_predictPositions.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_predictDensity(){
		// Stage predict density
		_pcisph_predictDensity.setArg( 0, _neighborMap );
		_pcisph_predictDensity.setArg( 1, _particleIndexBack );
		_pcisph_predictDensity.setArg( 2, PhysicsConstants.W_POLY_6_COEFFICIENT );
		_pcisph_predictDensity.setArg( 3, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		_pcisph_predictDensity.setArg( 4, PhysicsConstants.H );
		_pcisph_predictDensity.setArg( 5, PhysicsConstants.MASS );
		_pcisph_predictDensity.setArg( 6, PhysicsConstants.RHO0 );
		_pcisph_predictDensity.setArg( 7, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_predictDensity.setArg( 8, PhysicsConstants.STIFFNESS );
		_pcisph_predictDensity.setArg( 9, _sortedPosition );
		_pcisph_predictDensity.setArg(10, _pressure );
		_pcisph_predictDensity.setArg(11, _rho );
		_pcisph_predictDensity.setArg(12, PhysicsConstants.DELTA );
		_pcisph_predictDensity.setArg( 13, _particleCount );
		_pcisph_predictDensity.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_correctPressure(){
		// Stage correct pressure
		_pcisph_correctPressure.setArg( 0, _neighborMap );
		_pcisph_correctPressure.setArg( 1, _particleIndexBack );
		_pcisph_correctPressure.setArg( 2, PhysicsConstants.W_POLY_6_COEFFICIENT );
		_pcisph_correctPressure.setArg( 3, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		_pcisph_correctPressure.setArg( 4, PhysicsConstants.H );
		_pcisph_correctPressure.setArg( 5, PhysicsConstants.MASS );
		_pcisph_correctPressure.setArg( 6, PhysicsConstants.RHO0 );
		_pcisph_correctPressure.setArg( 7, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_correctPressure.setArg( 8, PhysicsConstants.STIFFNESS );
		_pcisph_correctPressure.setArg( 9, _sortedPosition );
		_pcisph_correctPressure.setArg(10, _pressure );
		_pcisph_correctPressure.setArg(11, _rho );
		_pcisph_correctPressure.setArg(12, PhysicsConstants.DELTA );
		_pcisph_correctPressure.setArg(13, _position );
		_pcisph_correctPressure.setArg(14, _particleIndex );
		_pcisph_correctPressure.setArg( 15, _particleCount );
		_pcisph_correctPressure.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_computePressureForceAcceleration(){
		// Stage ComputeAcceleration
		_pcisph_computePressureForceAcceleration.setArg( 0, _neighborMap );
		_pcisph_computePressureForceAcceleration.setArg( 1, _pressure );
		_pcisph_computePressureForceAcceleration.setArg( 2, _rho );
		_pcisph_computePressureForceAcceleration.setArg( 3, _sortedPosition );
		_pcisph_computePressureForceAcceleration.setArg( 4, _sortedVelocity );
		_pcisph_computePressureForceAcceleration.setArg( 5, _particleIndexBack );
		_pcisph_computePressureForceAcceleration.setArg( 6, PhysicsConstants.CFLLimit );
		_pcisph_computePressureForceAcceleration.setArg( 7, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT );
		_pcisph_computePressureForceAcceleration.setArg( 8, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		_pcisph_computePressureForceAcceleration.setArg( 9, PhysicsConstants.H );
		_pcisph_computePressureForceAcceleration.setArg( 10, PhysicsConstants.MASS );
		_pcisph_computePressureForceAcceleration.setArg( 11, PhysicsConstants.MU );
		_pcisph_computePressureForceAcceleration.setArg( 12, PhysicsConstants.SIMULATION_SCALE );
		_pcisph_computePressureForceAcceleration.setArg( 13, _acceleration );
		_pcisph_computePressureForceAcceleration.setArg( 14, PhysicsConstants.RHO0 );
		_pcisph_computePressureForceAcceleration.setArg( 15, _position );
		_pcisph_computePressureForceAcceleration.setArg( 16, _particleIndex );
		_pcisph_computePressureForceAcceleration.setArg( 17, _particleCount );
		_pcisph_computePressureForceAcceleration.enqueueNDRange(_queue, new int[] {_particleCount});
		
		return 0;
	}
	
	public int run_pcisph_integrate(){
		// Stage Integrate
		_pcisph_integrate.setArg( 0, _acceleration );
		_pcisph_integrate.setArg( 1, _sortedPosition );
		_pcisph_integrate.setArg( 2, _sortedVelocity );
		_pcisph_integrate.setArg( 3, _particleIndex );
		_pcisph_integrate.setArg( 4, _particleIndexBack );
		_pcisph_integrate.setArg( 5, PhysicsConstants.GRAVITY_X );
		_pcisph_integrate.setArg( 6, PhysicsConstants.GRAVITY_Y );
		_pcisph_integrate.setArg( 7, PhysicsConstants.GRAVITY_Z );
		_pcisph_integrate.setArg( 8, PhysicsConstants.SIMULATION_SCALE_INV );
		_pcisph_integrate.setArg( 9, PhysicsConstants.TIME_STEP );
		_pcisph_integrate.setArg( 10, SPHConstants.XMIN );
		_pcisph_integrate.setArg( 11, SPHConstants.XMAX );
		_pcisph_integrate.setArg( 12, SPHConstants.YMIN );
		_pcisph_integrate.setArg( 13, SPHConstants.YMAX );
		_pcisph_integrate.setArg( 14, SPHConstants.ZMIN );
		_pcisph_integrate.setArg( 15, SPHConstants.ZMAX );
		_pcisph_integrate.setArg( 16, PhysicsConstants.DAMPING );
		_pcisph_integrate.setArg( 17, _position );
		_pcisph_integrate.setArg( 18, _velocity );
		_pcisph_integrate.setArg( 19, _rho );
		_pcisph_integrate.setArg( 20, PhysicsConstants.R0 );
		_pcisph_integrate.setArg( 21, _neighborMap );
		_pcisph_integrate.setArg(22, _particleCount );
		_pcisph_integrate.enqueueNDRange(_queue, new int[] {_particleCount});
		
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
		logger.info("SPH clear buffer");
		runClearBuffers();
		logger.info("SPH hash particles");
		runHashParticles();
		logger.info("SPH sort");
		runSort();
		logger.info("SPH sort post pass");
		runSortPostPass();
		logger.info("SPH index");
		runIndexx();
		logger.info("SPH index post pass");
		runIndexPostPass();
		logger.info("SPH find neighbors");
		runFindNeighbors();
		
		// PORTING-NOTE: on the original C++ version preparing elastic matter stuff was done at this point
		// For this implementation we are moving elastic matter data in the configuration file for the model

		// TODO: PCISPH stuff goes here
		logger.info("PCI-SPH compute density");
		run_pcisph_computeDensity();
		logger.info("PCI-SPH compute forces and init pressure");
		run_pcisph_computeForcesAndInitPressure();
		
		if(_numOfElasticP > 0){
			logger.info("PCI-SPH compute elastic forces");
			run_pcisph_computeElasticForces();
		}
		
		logger.info("PCI-SPH predict/correct loop");
		// LOOP: 3 times or until "error" becomes less than 2%
		int iter = 0; int maxIterations = 3;
		do
		{
			run_pcisph_predictPositions();
			run_pcisph_predictDensity();
			run_pcisph_correctPressure();
			run_pcisph_computePressureForceAcceleration();
			
			iter++;
		} while ((iter<maxIterations));
		
		logger.info("PCI-SPH integrate");
		run_pcisph_integrate();
		
		logger.info("SPH position read");
		logger.info("Position element size: " + _position.getElementSize());
		logger.info("Position element count: " + _position.getElementCount());		
		_positionPtr = _position.read(_queue);
		
		logger.info("Position element size: " + _velocity.getElementSize());
		logger.info("Position element count: " + _velocity.getElementCount());
		_velocityPtr = _velocity.read(_queue);
		logger.info("SPH finish queue");
		_queue.finish();
		logger.info("SPH step done");
	}
	
	public void finishQueue() {
		_queue.finish();
	}
	
	public Float round(Float val, int roundingFactor){
		return (float) Math.round(val * roundingFactor) / roundingFactor;
	}
};