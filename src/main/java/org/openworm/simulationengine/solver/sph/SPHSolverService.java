package org.openworm.simulationengine.solver.sph;

import static java.lang.System.out;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.bridj.Pointer;
import org.openworm.simulationengine.core.constants.PhysicsConstants;
import org.openworm.simulationengine.core.model.IModel;
import org.openworm.simulationengine.core.model.MathUtils;
import org.openworm.simulationengine.core.simulation.ITimeConfiguration;
import org.openworm.simulationengine.core.solver.ISolver;
import org.openworm.simulationengine.model.sph.SPHConstants;
import org.openworm.simulationengine.model.sph.SPHModel;
import org.openworm.simulationengine.model.sph.SPHParticle;
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
	
	private CLContext context;
	public CLQueue queue;
	private CLProgram program;
	private CLDevice device;
	public CLBuffer<Float> acceleration;
	public CLBuffer<Integer> gridCellIndex;
	public CLBuffer<Integer> gridCellIndexFixedUp;
	public CLBuffer<Float> neighborMap;
	public CLBuffer<Integer> particleIndex;
	public CLBuffer<Float> position;
	public CLBuffer<Float> pressure;
	public CLBuffer<Float> rho;
	public CLBuffer<Float> rhoInv;
	public CLBuffer<Float> sortedPosition;
	public CLBuffer<Float> sortedVelocity;
	public CLBuffer<Float> velocity;
	
	public Pointer<Float> accelerationPtr;
	public Pointer<Integer> gridCellIndexPtr;
	public Pointer<Integer> gridCellIndexFixedUpPtr;
	public Pointer<Float> neighborMapPtr;
	public Pointer<Integer> particleIndexPtr;
	public Pointer<Float> positionPtr;
	public Pointer<Float> positionPtrbuff;
	public Pointer<Float> pressurePtr;
	public Pointer<Float> rhoPtr;
	public Pointer<Float> rhoInvPtr;
	public Pointer<Float> sortedPositionPtr;
	public Pointer<Float> sortedVelocityPtr;
	public Pointer<Float> velocityPtr;
	
	private CLKernel clearBuffers;
	private CLKernel computeAcceleration;
	private CLKernel computeDensityPressure;
	private CLKernel findNeighbors;
	private CLKernel hashParticles;
	private CLKernel indexPostPass;
	private CLKernel indexx;
	private CLKernel integrate;
	private CLKernel sortPostPass;
	
	SPHModel model;
	public int gridCellsX;
	public int gridCellsY;
	public int gridCellsZ;
	public int gridCellCount;
	
	public static Random randomGenerator = new Random();
	
	public SPHSolverService() throws Exception{
		this.init();
	}
		
	private void init() throws IOException  
	{
		gridCellsX = (int)( ( SPHConstants.XMAX - SPHConstants.XMIN ) / PhysicsConstants.H ) + 1;
		gridCellsY = (int)( ( SPHConstants.YMAX - SPHConstants.YMIN ) / PhysicsConstants.H ) + 1;
		gridCellsZ = (int)( ( SPHConstants.ZMAX - SPHConstants.ZMIN ) / PhysicsConstants.H ) + 1;
		gridCellCount = gridCellsX * gridCellsY * gridCellsZ;
		model = new SPHModelX(gridCellsX, gridCellsY, gridCellsZ);
		context = JavaCL.createBestContext(DeviceFeature.GPU);
		
		out.println("created "+ context);
		// an array with available devices
		CLDevice[] devices = context.getDevices();

		for(int i=0; i<devices.length; i++)
		{
			out.println("device - " + i + ": " + devices[i]);
		}	

		// have a look at the output and select a device
		device = devices[0];
		out.println("Version " + device.getOpenCLVersion());
		out.println("Version " + device.getDriverVersion());
		out.println("using "+ device);
		out.println("max workgroup size: " + device.getMaxWorkGroupSize());
		out.println("max workitems size: " + device.getMaxWorkItemSizes()[0]);
		
		// create command queue on selected device.
		queue = context.createDefaultQueue();//device.createCommandQueue();
		
		// load sources, create and build program
		String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/resource/sphFluidDemo.cl"));
		program = context.createProgram(src);//.build();
		
		/* input buffers declarations Pointer Alternative*/
		accelerationPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		gridCellIndexPtr = Pointer.allocateInts((gridCellCount + 1));
		gridCellIndexFixedUpPtr = Pointer.allocateInts((gridCellCount + 1));
		neighborMapPtr = Pointer.allocateFloats(SPHConstants.NK * 2);
		particleIndexPtr = Pointer.allocateInts(SPHConstants.PARTICLE_COUNT * 2);
		positionPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		pressurePtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 1);
		rhoPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 1);
		rhoInvPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 1);
		sortedPositionPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		sortedVelocityPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		velocityPtr = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		positionPtrbuff = Pointer.allocateFloats(SPHConstants.PARTICLE_COUNT * 4);
		/*end declaration*/
		
		/*kernels*/
		clearBuffers = program.createKernel("clearBuffers");
		computeAcceleration = program.createKernel("computeAcceleration");
		computeDensityPressure = program.createKernel("computeDensityPressure");
		findNeighbors = program.createKernel("findNeighbors");
		hashParticles = program.createKernel("hashParticles");
		indexPostPass = program.createKernel("indexPostPass");
		indexx = program.createKernel("indexx");
		integrate = program.createKernel("integrate");
		sortPostPass = program.createKernel("sortPostPass");
		
		/* particle creation */
		// TODO: move this out of here (unit tests)
		int index = 0;
		for(int i = 0;i<SPHConstants.PARTICLE_COUNT;i++){
			if(i != 0)
				index = index + 4;
			float r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			Vector3DX positionVector = new Vector3DX();
			Vector3DX velocityVector = new Vector3DX();
			positionVector.setX(MathUtils.scale(SPHConstants.XMIN, SPHConstants.XMAX/10 , r)); 
			r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			positionVector.setY(MathUtils.scale(SPHConstants.YMIN, SPHConstants.YMAX , r)); 
			r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			positionVector.setZ(MathUtils.scale(SPHConstants.ZMIN, SPHConstants.ZMAX , r));
			positionVector.setP(0f);
			positionPtr.set(index,positionVector.getX());
			positionPtr.set(index + 1,positionVector.getY());
			positionPtr.set(index + 2,positionVector.getZ());
			positionPtr.set(index + 3,positionVector.getP());
			r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			velocityVector.setX(MathUtils.scale(-1.0f, 1.0f, r));
			r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			velocityVector.setY(MathUtils.scale(-1.0f, 1.0f, r));
			r = ((float)randomGenerator.nextInt(PhysicsConstants.RAND_MAX) / (float)PhysicsConstants.RAND_MAX );//tr.rand.next();
			velocityVector.setZ(MathUtils.scale(-1.0f, 1.0f, r));
			velocityVector.setP(0f);
			velocityPtr.set(index,velocityVector.getX());
			velocityPtr.set(index + 1,velocityVector.getY());
			velocityPtr.set(index + 2,velocityVector.getZ());
			velocityPtr.set(index + 3,velocityVector.getP());
			SPHParticle particle = new SPHParticleX(positionVector, velocityVector, 1);
			model.getParticles().add(particle);
		}
		
		/*Alternativ buffer defining*/
		acceleration = context.createBuffer(Usage.InputOutput,accelerationPtr,false);
		gridCellIndex = context.createBuffer(Usage.InputOutput,gridCellIndexPtr, false);
		gridCellIndexFixedUp = context.createIntBuffer(Usage.InputOutput,gridCellIndexFixedUpPtr, false);
		neighborMap = context.createBuffer(Usage.InputOutput,neighborMapPtr, false);
		particleIndex = context.createBuffer(Usage.InputOutput,particleIndexPtr ,false);
		position = context.createBuffer(Usage.InputOutput, positionPtr, false);
		pressure = context.createBuffer(Usage.InputOutput,pressurePtr,false);
		rho = context.createBuffer(Usage.InputOutput,rhoPtr,false);
		rhoInv = context.createBuffer(Usage.InputOutput,rhoInvPtr,false);
		sortedPosition = context.createBuffer(Usage.InputOutput,sortedPositionPtr,false);
		sortedVelocity = context.createBuffer(Usage.InputOutput, sortedVelocityPtr,false);
		velocity = context.createBuffer(Usage.InputOutput,velocityPtr,false);
		queue.finish();
	}
	
	public List<List<IModel>> solve(List<IModel> models, ITimeConfiguration timeConfiguration)
	{
		// TODO: extend this to use time configuration to do multiple steps in one go
		
		// TODO: 1. populate buffers from list of models
		// TODO: 2. call this.step();
		// TODO: 3. retrieve values from buffers and populate returned models
		return null;
	}
	
	public void cleanContext(){
		this.context.release();
	}
	
	public int _runClearBuffers(){
		clearBuffers.setArg(0, neighborMap);
		clearBuffers.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runComputeAcceleration(){
		computeAcceleration.setArg( 0, neighborMap );
		computeAcceleration.setArg( 1, pressure );
		computeAcceleration.setArg( 2, rho );
		computeAcceleration.setArg( 3, rhoInv );
		computeAcceleration.setArg( 4, sortedPosition );
		computeAcceleration.setArg( 5, sortedVelocity );
		computeAcceleration.setArg( 6, PhysicsConstants.CFLLimit );
		computeAcceleration.setArg( 7, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT );
		computeAcceleration.setArg( 8, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT );
		computeAcceleration.setArg( 9, PhysicsConstants.H );
		computeAcceleration.setArg( 10, PhysicsConstants.MASS );
		computeAcceleration.setArg( 11, PhysicsConstants.MU );
		computeAcceleration.setArg( 12, PhysicsConstants.SIMULATION_SCALE );
		computeAcceleration.setArg( 13, acceleration );
		computeAcceleration.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runComputeDensityPressure(){
		computeDensityPressure.setArg( 0, neighborMap );
		computeDensityPressure.setArg( 1, PhysicsConstants.W_POLY_6_COEFFICIENT );
		computeDensityPressure.setArg( 2, PhysicsConstants.H );
		computeDensityPressure.setArg( 3, PhysicsConstants.MASS );
		computeDensityPressure.setArg( 4, PhysicsConstants.RHO0 );
		computeDensityPressure.setArg( 5, PhysicsConstants.SIMULATION_SCALE );
		computeDensityPressure.setArg( 6, PhysicsConstants.STIFFNESS );
		computeDensityPressure.setArg( 7, pressure );
		computeDensityPressure.setArg( 8, rho );
		computeDensityPressure.setArg( 9, rhoInv );
		computeDensityPressure.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runFindNeighbors(){
		findNeighbors.setArg( 0, gridCellIndexFixedUp );
		findNeighbors.setArg( 1, sortedPosition );
		gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
		findNeighbors.setArg( 2, gridCellCount );
		findNeighbors.setArg( 3, gridCellsX );
		findNeighbors.setArg( 4, gridCellsY );
		findNeighbors.setArg( 5, gridCellsZ );
		findNeighbors.setArg( 6, PhysicsConstants.H );
		findNeighbors.setArg( 7, PhysicsConstants.HASH_GRID_CELL_SIZE );
		findNeighbors.setArg( 8, PhysicsConstants.HASH_GRID_CELL_SIZE_INV );
		findNeighbors.setArg( 9, PhysicsConstants.SIMULATION_SCALE );
		findNeighbors.setArg( 10, SPHConstants.XMIN_F );
		findNeighbors.setArg( 11, SPHConstants.YMIN_F );
		findNeighbors.setArg( 12, SPHConstants.ZMIN_F );
		findNeighbors.setArg( 13, neighborMap );
		findNeighbors.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runHashParticles(){
		// Stage HashParticles
		hashParticles.setArg( 0, position );
		hashParticles.setArg( 1, gridCellsX );
		hashParticles.setArg( 2, gridCellsY );
		hashParticles.setArg( 3, gridCellsZ );
		hashParticles.setArg( 4, PhysicsConstants.HASH_GRID_CELL_SIZE_INV );
		hashParticles.setArg( 5, SPHConstants.XMIN_F );
		hashParticles.setArg( 6, SPHConstants.YMIN_F );
		hashParticles.setArg( 7, SPHConstants.ZMIN_F );
		hashParticles.setArg( 8, particleIndex );
		hashParticles.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runIndexPostPass(){
		indexPostPass.setArg( 0, gridCellIndex );
		gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
		indexPostPass.setArg( 1, gridCellCount );
		indexPostPass.setArg( 2, gridCellIndexFixedUp );
		int gridCellCountRoundedUp = ((( gridCellCount - 1 ) / 256 ) + 1 ) * 256;
		indexPostPass.enqueueNDRange(queue, new int[] {gridCellCountRoundedUp});
		return 0;
	}
	
	public int _runIndexx(){
		// Stage Indexx
		indexx.setArg( 0, particleIndex );
		gridCellCount = ((gridCellsX) * (gridCellsY)) * (gridCellsZ);
		indexx.setArg( 1, gridCellCount );
		indexx.setArg( 2, gridCellIndex );
		int gridCellCountRoundedUp = ((( gridCellCount - 1 ) / 256 ) + 1 ) * 256;
		indexx.enqueueNDRange(queue, new int[] {gridCellCountRoundedUp});
		return 0;
	}
	
	public int _runIntegrate(){
		// Stage Integrate
		integrate.setArg( 0, acceleration );
		integrate.setArg( 1, sortedPosition );
		integrate.setArg( 2, sortedVelocity );
		integrate.setArg( 3, PhysicsConstants.GRAVITY_X );
		integrate.setArg( 4, PhysicsConstants.GRAVITY_Y );
		integrate.setArg( 5, PhysicsConstants.GRAVITY_Z );
		integrate.setArg( 6, PhysicsConstants.SIMULATION_SCALE_INV );
		integrate.setArg( 7, PhysicsConstants.TIME_STEP );
		integrate.setArg( 8, SPHConstants.XMIN_F );
		integrate.setArg( 9, SPHConstants.XMAX_F );
		integrate.setArg( 10, SPHConstants.YMIN_F );
		integrate.setArg( 11, SPHConstants.YMAX_F );
		integrate.setArg( 12, SPHConstants.ZMIN_F );
		integrate.setArg( 13, SPHConstants.ZMAX_F );
		integrate.setArg( 14, PhysicsConstants.DAMPING );
		integrate.setArg( 15, position );
		integrate.setArg( 16, velocity );
		integrate.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runSortPostPass(){
		// Stage SortPostPass
		sortPostPass.setArg( 0, particleIndex );
		sortPostPass.setArg( 1, position );
		sortPostPass.setArg( 2, velocity );
		sortPostPass.setArg( 3, sortedPosition );
		sortPostPass.setArg( 4, sortedVelocity );
		sortPostPass.enqueueNDRange(queue, new int[] {SPHConstants.PARTICLE_COUNT});
		return 0;
	}
	
	public int _runSort(){
		//this version work with qsort
		int index = 0;
		List<int[]> _particleIndex = new ArrayList<int[]>();
		Pointer<Integer> particleInd = particleIndex.read(queue);
		queue.finish();
		for(int i = 0; i < SPHConstants.PARTICLE_COUNT * 2;i+=2){
			int[] element = {particleInd.get(i), particleInd.get(i+1)};
			_particleIndex.add(element);
		}
		Collections.sort(_particleIndex, new MyCompare());
		for(int i = 0; i< _particleIndex.size();i++){
			for(int j=0;j<2;j++){
				particleInd.set(index,_particleIndex.get(i)[j]);
				index++;
			}
		}
		particleIndex.write(queue, particleInd, false);
		queue.finish();
		return 0;
	}
	
	class MyCompare implements Comparator<int[]>{
		public int compare(int[] o1, int[] o2) {
			if( o1[0] < o2[0] ) return -1;
			if( o1[0] > o2[0]) return +1;
			return 0;
		}
	}
	
	void advanceInTime(){
		step();
		for (int i = 0;i<model.getParticles().size();i++) {
			model.getParticles().get(i).getPositionVector().setX(positionPtr.get(4*i));
			model.getParticles().get(i).getPositionVector().setY(positionPtr.get(4*i+1));
			model.getParticles().get(i).getPositionVector().setZ(positionPtr.get(4*i+2));
			model.getParticles().get(i).getPositionVector().setP(positionPtr.get(4*i+3));
		}
	}
	
	private void step(){
		_runClearBuffers();
		_runHashParticles();
		_runSort();
		_runSortPostPass();
		_runIndexx();
		_runIndexPostPass();
		_runFindNeighbors();
		_runComputeDensityPressure();
		_runComputeAcceleration();
		_runIntegrate();
		positionPtr = position.read(queue);
		queue.finish();
	}
	
	public void finishQueue() {
		queue.finish();
	}
}
;