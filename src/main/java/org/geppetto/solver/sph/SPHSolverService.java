package org.geppetto.solver.sph;

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
import org.geppetto.core.common.GeppettoInitializationException;
import org.geppetto.core.constants.PhysicsConstants;
import org.geppetto.core.model.IModel;
import org.geppetto.core.model.StateInstancePath;
import org.geppetto.core.model.StateSet;
import org.geppetto.core.model.values.ValuesFactory;
import org.geppetto.core.simulation.IRunConfiguration;
import org.geppetto.core.solver.ISolver;
import org.geppetto.model.sph.Connection;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.Vector3DX;
import org.springframework.stereotype.Service;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem;
import com.nativelibs4java.opencl.CLPlatform.DeviceFeature;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.util.IOUtils;

@Service
public class SPHSolverService implements ISolver
{

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

	public float _xMax;
	public float _xMin;
	public float _yMax;
	public float _yMin;
	public float _zMax;
	public float _zMin;
	public float _muscle_activation_signal = 0.0f;

	public int _gridCellsX;
	public int _gridCellsY;
	public int _gridCellsZ;
	public int _gridCellCount;
	public int _particleCount;
	public int _numOfLiquidP;
	public int _numOfElasticP;
	public int _numOfBoundaryP;

	private SPHModelX _model;

	public static Random RandomGenerator = new Random();

	public SPHSolverService(String hardwareProfile) throws Exception
	{
		this.onceOffInit(hardwareProfile);
	}
	
	public SPHSolverService() throws Exception
	{
		this(SPHConstants.CPU_PROFILE);
	}

	private void onceOffInit(String hwProfile) throws IOException
	{
		// TODO: check if the selected profile is actually available
		DeviceFeature feature = (hwProfile == SPHConstants.CPU_PROFILE) ? DeviceFeature.CPU : DeviceFeature.GPU;
		
		_context = JavaCL.createBestContext(feature);

		out.println("created " + _context);
		// an array with available devices
		CLDevice[] devices = _context.getDevices();

		for(int i = 0; i < devices.length; i++)
		{
			out.println("device - " + i + ": " + devices[i]);
		}

		// have a look at the output and select a device
		_device = devices[0];
		out.println("Version " + _device.getOpenCLVersion());
		out.println("Version " + _device.getDriverVersion());
		out.println("using " + _device);
		out.println("max workgroup size: " + _device.getMaxWorkGroupSize());
		out.println("max workitems size: " + _device.getMaxWorkItemSizes()[0]);

		// create command queue on selected device.
		_queue = _context.createDefaultQueue();// device.createCommandQueue();

		// load sources, create and build program
		String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/resource/sphFluid.cl"));
		_program = _context.createProgram(src);

		// kernels
		_clearBuffers = _program.createKernel("clearBuffers");
		_findNeighbors = _program.createKernel("findNeighbors");
		_hashParticles = _program.createKernel("hashParticles");
		// PORTING-NOTE: indexPostPass is gone from kernels in the latest version,logic moved from opencl to host code (this class)
		_indexx = _program.createKernel("indexx");
		_sortPostPass = _program.createKernel("sortPostPass");

		// PCI-SPH specific
		_pcisph_computeForcesAndInitPressure = _program.createKernel("pcisph_computeForcesAndInitPressure");
		_pcisph_integrate = _program.createKernel("pcisph_integrate");
		_pcisph_predictPositions = _program.createKernel("pcisph_predictPositions");
		_pcisph_predictDensity = _program.createKernel("pcisph_predictDensity");
		_pcisph_correctPressure = _program.createKernel("pcisph_correctPressure");
		_pcisph_computePressureForceAcceleration = _program.createKernel("pcisph_computePressureForceAcceleration");
		_pcisph_computeDensity = _program.createKernel("pcisph_computeDensity");
		_pcisph_computeElasticForces = _program.createKernel("pcisph_computeElasticForces");
	}

	private void allocateBuffers(){		
		// allocate native device memory for all buffers
		_acceleration = _context.createFloatBuffer(CLMem.Usage.InputOutput, _particleCount * 4 * 2);
		_gridCellIndex = _context.createIntBuffer(CLMem.Usage.InputOutput, _gridCellCount + 1);
		_gridCellIndexFixedUp = _context.createIntBuffer(CLMem.Usage.Input, _gridCellCount + 1);
		_neighborMap = _context.createFloatBuffer(CLMem.Usage.Input, _particleCount * SPHConstants.NEIGHBOR_COUNT * 2);
		_particleIndex = _context.createIntBuffer(CLMem.Usage.InputOutput, _particleCount * 2);
		_particleIndexBack = _context.createIntBuffer(CLMem.Usage.Input, _particleCount);
		_position = _context.createFloatBuffer(CLMem.Usage.InputOutput, _particleCount * 4);
		_pressure = _context.createFloatBuffer(CLMem.Usage.Input, _particleCount);
		_rho = _context.createFloatBuffer(CLMem.Usage.Input, _particleCount * 2);
		_sortedPosition = _context.createFloatBuffer(CLMem.Usage.Input, _particleCount * 4 * 2);
		_sortedVelocity = _context.createFloatBuffer(CLMem.Usage.Input, _particleCount * 4);
		_velocity = _context.createFloatBuffer(CLMem.Usage.InputOutput, _particleCount * 4);
	}

	public void setBuffersFromModel()
	{
		// set dimensions
		_xMax = _model.getXMax();
		_xMin = _model.getXMin();
		_yMax = _model.getYMax();
		_yMin = _model.getYMin();
		_zMax = _model.getZMax();
		_zMin = _model.getZMin();

		_particleCount = _model.getNumberOfParticles();
		_numOfElasticP = 0;
		_numOfLiquidP = 0;
		_numOfBoundaryP = 0;

		_gridCellsX = (int) ((_model.getXMax() - _model.getXMin()) / PhysicsConstants.H) + 1;
		_gridCellsY = (int) ((_model.getYMax() - _model.getYMin()) / PhysicsConstants.H) + 1;
		_gridCellsZ = (int) ((_model.getZMax() - _model.getZMin()) / PhysicsConstants.H) + 1;

		// set grid dimensions
		_gridCellCount = _gridCellsX * _gridCellsY * _gridCellsZ;

		// allocate buffers - requires global dimensions of the grid
		this.allocateBuffers();

		int index = 0;

		for(int i = 0; i < _particleCount; i++)
		{
			if(i != 0)
			{
				index = index + 4;
			}

			Vector3DX positionVector = (Vector3DX) _model.getParticles().get(i).getPositionVector();
			Vector3DX velocityVector = (Vector3DX) _model.getParticles().get(i).getVelocityVector();
			
			// map for writing
			_positionPtr = _position.map(_queue, CLMem.MapFlags.Write);
			_velocityPtr = _velocity.map(_queue, CLMem.MapFlags.Write);
			
			// buffer population
			_positionPtr.set(index, positionVector.getX());
			_positionPtr.set(index + 1, positionVector.getY());
			_positionPtr.set(index + 2, positionVector.getZ());
			_positionPtr.set(index + 3, positionVector.getP());
			_velocityPtr.set(index, velocityVector.getX());
			_velocityPtr.set(index + 1, velocityVector.getY());
			_velocityPtr.set(index + 2, velocityVector.getZ());
			_velocityPtr.set(index + 3, velocityVector.getP());
			
			// unmap after writing
			_position.unmap(_queue, _positionPtr);
			_velocity.unmap(_queue, _velocityPtr);

			// particle counts
			if(positionVector.getP() == SPHConstants.BOUNDARY_TYPE)
			{
				_numOfBoundaryP++;
			}
			else if(positionVector.getP() == SPHConstants.ELASTIC_TYPE)
			{
				_numOfElasticP++;
			}
			else if(positionVector.getP() == SPHConstants.LIQUID_TYPE)
			{
				_numOfLiquidP++;
			}
		}

		// populate elastic connection buffers if we have any
		if(_numOfElasticP > 0 && _model.getConnections().size() > 0)
		{
			// init elastic connections buffers
			// TODO: move this back with the other buffers init stuff
			_elasticConnectionsData = _context.createFloatBuffer(CLMem.Usage.InputOutput, _numOfElasticP * SPHConstants.NEIGHBOR_COUNT * 4);
			_elasticConnectionsDataPtr = _elasticConnectionsData.map(_queue, CLMem.MapFlags.Write);
				
			int connIndex = 0;
			for(Connection conn : _model.getConnections())
			{
				_elasticConnectionsDataPtr.set(connIndex, conn.getP1());
				_elasticConnectionsDataPtr.set(connIndex + 1, conn.getP2());
				_elasticConnectionsDataPtr.set(connIndex + 2, conn.getDistance());
				_elasticConnectionsDataPtr.set(connIndex + 3, 0f); // padding
				connIndex++;
			}
			
			// we copied the stuff down to the device and we won't touch it again so we can unmap
			_elasticConnectionsData.unmap(_queue, _elasticConnectionsDataPtr);
		}

		// check that counts are fine
		if(_particleCount != (_numOfBoundaryP + _numOfElasticP + _numOfLiquidP))
		{
			throw new IllegalArgumentException("SPHSolverService:setModels - particle counts do not add up");
		}
	}

	public void cleanContext()
	{
		this._context.release();
	}

	public int runClearBuffers()
	{
		_clearBuffers.setArg(0, _neighborMap);
		_clearBuffers.setArg(1, _particleCount);
		_clearBuffers.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	public int runFindNeighbors()
	{
		_findNeighbors.setArg(0, _gridCellIndexFixedUp);
		_findNeighbors.setArg(1, _sortedPosition);
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_findNeighbors.setArg(2, _gridCellCount);
		_findNeighbors.setArg(3, _gridCellsX);
		_findNeighbors.setArg(4, _gridCellsY);
		_findNeighbors.setArg(5, _gridCellsZ);
		_findNeighbors.setArg(6, PhysicsConstants.H);
		_findNeighbors.setArg(7, PhysicsConstants.HASH_GRID_CELL_SIZE);
		_findNeighbors.setArg(8, PhysicsConstants.HASH_GRID_CELL_SIZE_INV);
		_findNeighbors.setArg(9, PhysicsConstants.SIMULATION_SCALE);
		_findNeighbors.setArg(10, _xMin);
		_findNeighbors.setArg(11, _yMin);
		_findNeighbors.setArg(12, _zMin);
		_findNeighbors.setArg(13, _neighborMap);
		_findNeighbors.setArg(14, _particleCount);
		_findNeighbors.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	public CLEvent runHashParticles(){
		// Stage HashParticles
		_hashParticles.setArg( 0, _position );
		_hashParticles.setArg( 1, _gridCellsX );
		_hashParticles.setArg( 2, _gridCellsY );
		_hashParticles.setArg( 3, _gridCellsZ );
		_hashParticles.setArg( 4, PhysicsConstants.HASH_GRID_CELL_SIZE_INV );
		_hashParticles.setArg( 5, _xMin );
		_hashParticles.setArg( 6, _yMin );
		_hashParticles.setArg( 7, _zMin );
		_hashParticles.setArg( 8, _particleIndex );
		_hashParticles.setArg( 9, _particleCount );
		CLEvent event = _hashParticles.enqueueNDRange(_queue, new int[] {getParticleCountRoundedUp()});
		
		return event;
	}

	public int runIndexPostPass(){		
		// get values out of buffer 
		_gridCellIndexPtr = _gridCellIndex.map(_queue, CLMem.MapFlags.Read);
		int[] gridNextNonEmptyCellBuffer = _gridCellIndexPtr.getInts();
		_gridCellIndex.unmap(_queue, _gridCellIndexPtr);
		
		int recentNonEmptyCell = _gridCellCount;
		for(int i= _gridCellCount; i>=0; i--)
		{
			if(gridNextNonEmptyCellBuffer[i] == SPHConstants.NO_CELL_ID) {
				gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell; 
			}
			else {
				recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
			}
		}
		
		// put results back
		_gridCellIndexFixedUpPtr = _gridCellIndexFixedUp.map(_queue, CLMem.MapFlags.Write);
		_gridCellIndexFixedUpPtr.setInts(gridNextNonEmptyCellBuffer);
		_gridCellIndexFixedUp.unmap(_queue, _gridCellIndexFixedUpPtr);
		
		return 0;
	}

	public CLEvent runIndexx(){
		// Stage Indexx
		_indexx.setArg( 0, _particleIndex );
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_indexx.setArg( 1, _gridCellCount );
		_indexx.setArg( 2, _gridCellIndex );
		_indexx.setArg( 3, _particleCount );
		int gridCellCountRoundedUp = ((( _gridCellCount - 1 ) / 256 ) + 1 ) * 256;
		CLEvent event = _indexx.enqueueNDRange(_queue, new int[] {gridCellCountRoundedUp});
		
		return event;
	}

	public int runSortPostPass()
	{
		// Stage SortPostPass
		_sortPostPass.setArg(0, _particleIndex);
		_sortPostPass.setArg(1, _particleIndexBack);
		_sortPostPass.setArg(2, _position);
		_sortPostPass.setArg(3, _velocity);
		_sortPostPass.setArg(4, _sortedPosition);
		_sortPostPass.setArg(5, _sortedVelocity);
		_sortPostPass.setArg(6, _particleCount);
		_sortPostPass.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	public int run_pcisph_computeDensity()
	{
		// Stage ComputeDensityPressure
		_pcisph_computeDensity.setArg(0, _neighborMap);
		_pcisph_computeDensity.setArg(1, PhysicsConstants.W_POLY_6_COEFFICIENT);
		_pcisph_computeDensity.setArg(2, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_computeDensity.setArg(3, PhysicsConstants.H);
		_pcisph_computeDensity.setArg(4, PhysicsConstants.MASS);
		_pcisph_computeDensity.setArg(5, PhysicsConstants.RHO0);
		_pcisph_computeDensity.setArg(6, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_computeDensity.setArg(7, PhysicsConstants.STIFFNESS);
		_pcisph_computeDensity.setArg(8, _sortedPosition);
		_pcisph_computeDensity.setArg(9, _pressure);
		_pcisph_computeDensity.setArg(10, _rho);
		_pcisph_computeDensity.setArg(11, _particleIndexBack);
		_pcisph_computeDensity.setArg(12, PhysicsConstants.DELTA); // calculated from constants
		_pcisph_computeDensity.setArg(13, _particleCount);
		_pcisph_computeDensity.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public int run_pcisph_computeForcesAndInitPressure()
	{
		_pcisph_computeForcesAndInitPressure.setArg(0, _neighborMap);
		_pcisph_computeForcesAndInitPressure.setArg(1, _rho);
		_pcisph_computeForcesAndInitPressure.setArg(2, _pressure);
		_pcisph_computeForcesAndInitPressure.setArg(3, _sortedPosition);
		_pcisph_computeForcesAndInitPressure.setArg(4, _sortedVelocity);
		_pcisph_computeForcesAndInitPressure.setArg(5, _acceleration);
		_pcisph_computeForcesAndInitPressure.setArg(6, _particleIndexBack);
		_pcisph_computeForcesAndInitPressure.setArg(7, PhysicsConstants.W_POLY_6_COEFFICIENT);
		_pcisph_computeForcesAndInitPressure.setArg(8, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT);
		_pcisph_computeForcesAndInitPressure.setArg(9, PhysicsConstants.H);
		_pcisph_computeForcesAndInitPressure.setArg(10, PhysicsConstants.MASS);
		_pcisph_computeForcesAndInitPressure.setArg(11, PhysicsConstants.MU);
		_pcisph_computeForcesAndInitPressure.setArg(12, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_computeForcesAndInitPressure.setArg(13, PhysicsConstants.GRAVITY_X);
		_pcisph_computeForcesAndInitPressure.setArg(14, PhysicsConstants.GRAVITY_Y);
		_pcisph_computeForcesAndInitPressure.setArg(15, PhysicsConstants.GRAVITY_Z);
		_pcisph_computeForcesAndInitPressure.setArg(16, _position);
		_pcisph_computeForcesAndInitPressure.setArg(17, _particleIndex);
		_pcisph_computeForcesAndInitPressure.setArg(18, _particleCount);
		_pcisph_computeForcesAndInitPressure.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public int run_pcisph_computeElasticForces()
	{
		_pcisph_computeElasticForces.setArg(0, _neighborMap);
		_pcisph_computeElasticForces.setArg(1, _sortedPosition);
		_pcisph_computeElasticForces.setArg(2, _sortedVelocity);
		_pcisph_computeElasticForces.setArg(3, _acceleration);
		_pcisph_computeElasticForces.setArg(4, _particleIndexBack);
		_pcisph_computeElasticForces.setArg(5, _velocity);
		_pcisph_computeElasticForces.setArg(6, PhysicsConstants.H);
		_pcisph_computeElasticForces.setArg(7, PhysicsConstants.MASS);
		_pcisph_computeElasticForces.setArg(8, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_computeElasticForces.setArg(9, _numOfElasticP);
		_pcisph_computeElasticForces.setArg(10, _elasticConnectionsData);
		_pcisph_computeElasticForces.setArg(11, _numOfBoundaryP);
		_pcisph_computeElasticForces.setArg(12, _muscle_activation_signal);
		_pcisph_computeElasticForces.setArg(13, _particleCount);

		int numOfElasticPRoundedUp = (((_numOfElasticP - 1) / 256) + 1) * 256;

		_pcisph_computeElasticForces.enqueueNDRange(_queue, new int[] { numOfElasticPRoundedUp });

		return 0;
	}

	public int run_pcisph_predictPositions()
	{
		_pcisph_predictPositions.setArg(0, _acceleration);
		_pcisph_predictPositions.setArg(1, _sortedPosition);
		_pcisph_predictPositions.setArg(2, _sortedVelocity);
		_pcisph_predictPositions.setArg(3, _particleIndex);
		_pcisph_predictPositions.setArg(4, _particleIndexBack);
		_pcisph_predictPositions.setArg(5, PhysicsConstants.GRAVITY_X);
		_pcisph_predictPositions.setArg(6, PhysicsConstants.GRAVITY_Y);
		_pcisph_predictPositions.setArg(7, PhysicsConstants.GRAVITY_Z);
		_pcisph_predictPositions.setArg(8, PhysicsConstants.SIMULATION_SCALE_INV);
		_pcisph_predictPositions.setArg(9, PhysicsConstants.TIME_STEP);
		_pcisph_predictPositions.setArg(10, _xMin);
		_pcisph_predictPositions.setArg(11, _xMax);
		_pcisph_predictPositions.setArg(12, _yMin);
		_pcisph_predictPositions.setArg(13, _yMax);
		_pcisph_predictPositions.setArg(14, _zMin);
		_pcisph_predictPositions.setArg(15, _zMax);
		_pcisph_predictPositions.setArg(16, PhysicsConstants.DAMPING);
		_pcisph_predictPositions.setArg(17, _position);
		_pcisph_predictPositions.setArg(18, _velocity);
		_pcisph_predictPositions.setArg(19, PhysicsConstants.R0);
		_pcisph_predictPositions.setArg(20, _neighborMap);
		_pcisph_predictPositions.setArg(21, _particleCount);
		_pcisph_predictPositions.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public int run_pcisph_predictDensity()
	{
		// Stage predict density
		_pcisph_predictDensity.setArg(0, _neighborMap);
		_pcisph_predictDensity.setArg(1, _particleIndexBack);
		_pcisph_predictDensity.setArg(2, PhysicsConstants.W_POLY_6_COEFFICIENT);
		_pcisph_predictDensity.setArg(3, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_predictDensity.setArg(4, PhysicsConstants.H);
		_pcisph_predictDensity.setArg(5, PhysicsConstants.MASS);
		_pcisph_predictDensity.setArg(6, PhysicsConstants.RHO0);
		_pcisph_predictDensity.setArg(7, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_predictDensity.setArg(8, PhysicsConstants.STIFFNESS);
		_pcisph_predictDensity.setArg(9, _sortedPosition);
		_pcisph_predictDensity.setArg(10, _pressure);
		_pcisph_predictDensity.setArg(11, _rho);
		_pcisph_predictDensity.setArg(12, PhysicsConstants.DELTA);
		_pcisph_predictDensity.setArg(13, _particleCount);
		_pcisph_predictDensity.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public int run_pcisph_correctPressure()
	{
		// Stage correct pressure
		_pcisph_correctPressure.setArg(0, _neighborMap);
		_pcisph_correctPressure.setArg(1, _particleIndexBack);
		_pcisph_correctPressure.setArg(2, PhysicsConstants.W_POLY_6_COEFFICIENT);
		_pcisph_correctPressure.setArg(3, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_correctPressure.setArg(4, PhysicsConstants.H);
		_pcisph_correctPressure.setArg(5, PhysicsConstants.MASS);
		_pcisph_correctPressure.setArg(6, PhysicsConstants.RHO0);
		_pcisph_correctPressure.setArg(7, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_correctPressure.setArg(8, PhysicsConstants.STIFFNESS);
		_pcisph_correctPressure.setArg(9, _sortedPosition);
		_pcisph_correctPressure.setArg(10, _pressure);
		_pcisph_correctPressure.setArg(11, _rho);
		_pcisph_correctPressure.setArg(12, PhysicsConstants.DELTA);
		_pcisph_correctPressure.setArg(13, _position);
		_pcisph_correctPressure.setArg(14, _particleIndex);
		_pcisph_correctPressure.setArg(15, _particleCount);
		_pcisph_correctPressure.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public int run_pcisph_computePressureForceAcceleration()
	{
		// Stage ComputeAcceleration
		_pcisph_computePressureForceAcceleration.setArg(0, _neighborMap);
		_pcisph_computePressureForceAcceleration.setArg(1, _pressure);
		_pcisph_computePressureForceAcceleration.setArg(2, _rho);
		_pcisph_computePressureForceAcceleration.setArg(3, _sortedPosition);
		_pcisph_computePressureForceAcceleration.setArg(4, _sortedVelocity);
		_pcisph_computePressureForceAcceleration.setArg(5, _particleIndexBack);
		_pcisph_computePressureForceAcceleration.setArg(6, PhysicsConstants.CFLLimit);
		_pcisph_computePressureForceAcceleration.setArg(7, PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT);
		_pcisph_computePressureForceAcceleration.setArg(8, PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_computePressureForceAcceleration.setArg(9, PhysicsConstants.H);
		_pcisph_computePressureForceAcceleration.setArg(10, PhysicsConstants.MASS);
		_pcisph_computePressureForceAcceleration.setArg(11, PhysicsConstants.MU);
		_pcisph_computePressureForceAcceleration.setArg(12, PhysicsConstants.SIMULATION_SCALE);
		_pcisph_computePressureForceAcceleration.setArg(13, _acceleration);
		_pcisph_computePressureForceAcceleration.setArg(14, PhysicsConstants.RHO0);
		_pcisph_computePressureForceAcceleration.setArg(15, _position);
		_pcisph_computePressureForceAcceleration.setArg(16, _particleIndex);
		_pcisph_computePressureForceAcceleration.setArg(17, _particleCount);
		_pcisph_computePressureForceAcceleration.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	public CLEvent run_pcisph_integrate(){
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
		_pcisph_integrate.setArg( 10, _xMin );
		_pcisph_integrate.setArg( 11, _xMax );
		_pcisph_integrate.setArg( 12, _yMin );
		_pcisph_integrate.setArg( 13, _yMax );
		_pcisph_integrate.setArg( 14, _zMin );
		_pcisph_integrate.setArg( 15, _zMax );
		_pcisph_integrate.setArg( 16, PhysicsConstants.DAMPING );
		_pcisph_integrate.setArg( 17, _position );
		_pcisph_integrate.setArg( 18, _velocity );
		_pcisph_integrate.setArg( 19, _rho );
		_pcisph_integrate.setArg( 20, PhysicsConstants.R0 );
		_pcisph_integrate.setArg( 21, _neighborMap );
		_pcisph_integrate.setArg(22, _particleCount );
		CLEvent event = _pcisph_integrate.enqueueNDRange(_queue, new int[] {getParticleCountRoundedUp()});
		
		return event;
	}

	public int runSort(){
		//this version work with qsort
		int index = 0;
		List<int[]> particleIndex = new ArrayList<int[]>();
		
		// get values out of buffer
		_particleIndexPtr = _particleIndex.map(_queue, CLMem.MapFlags.ReadWrite);
		int[] particleInd = _particleIndexPtr.getInts();
		
		for(int i = 0; i < _particleCount * 2;i+=2){
			int[] element = {particleInd[i], particleInd[i+1]};
			particleIndex.add(element);
		}
		Collections.sort(particleIndex, new MyCompare());
		for(int i = 0; i< particleIndex.size();i++){
			for(int j=0;j<2;j++){
				particleInd[index] = particleIndex.get(i)[j];
				index++;
			}
		}
		
		// put results back
		_particleIndexPtr.setInts(particleInd);
		_particleIndex.unmap(_queue, _particleIndexPtr);
		
		return 0;
	}

	class MyCompare implements Comparator<int[]>
	{
		public int compare(int[] o1, int[] o2)
		{
			if(o1[0] < o2[0]) return -1;
			if(o1[0] > o2[0]) return +1;
			return 0;
		}
	}

	private void step(){
		long endStep=0;
		long startStep=System.currentTimeMillis();
		long end=0;
		long start=System.currentTimeMillis();
		
		logger.info("SPH clear buffer");
		runClearBuffers();
		end=System.currentTimeMillis();
		logger.info("SPH clear buffer end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("SPH hash particles");
		CLEvent hashParticles = runHashParticles();
		end=System.currentTimeMillis();
		logger.info("SPH hash particles end, took "+ (end-start) +"ms");
		start=end;
		
		// host needs to wait as the next operation requires values from buffers
		hashParticles.waitFor();
		
		logger.info("SPH sort");
		runSort();
		end=System.currentTimeMillis();
		logger.info("SPH sort end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("SPH sort post pass");
		runSortPostPass();
		end=System.currentTimeMillis();
		logger.info("SPH sort post pass end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("SPH index");
		CLEvent runIndexx = runIndexx();
		end=System.currentTimeMillis();
		logger.info("SPH index end, took "+ (end-start) +"ms");
		start=end;
		
		// host needs to wait as the next operation requires values from buffers
		runIndexx.waitFor();
		
		logger.info("SPH index post pass");
		runIndexPostPass();
		end=System.currentTimeMillis();
		logger.info("SPH index post pass end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("SPH find neighbors");
		runFindNeighbors();
		end=System.currentTimeMillis();
		logger.info("SPH find neighbors end, took "+ (end-start) +"ms");
		start=end;
		
		// PCISPH stuff starts here
		logger.info("PCI-SPH compute density");
		run_pcisph_computeDensity();
		end=System.currentTimeMillis();
		logger.info("PCI-SPH compute density end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("PCI-SPH compute forces and init pressure");
		run_pcisph_computeForcesAndInitPressure();
		end=System.currentTimeMillis();
		logger.info("PCI-SPH compute forces and init pressure end, took "+ (end-start) +"ms");
		start=end;
		
		// Do elastic stuff only if we have elastic particles
		if(_numOfElasticP > 0){
			logger.info("PCI-SPH compute elastic forces");
			run_pcisph_computeElasticForces();
			end=System.currentTimeMillis();
			logger.info("PCI-SPH compute elastic forces end, took "+ (end-start) +"ms");
			start=end;
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
		end=System.currentTimeMillis();
		logger.info("PCI-SPH predict/correct loop end, took "+ (end-start) +"ms");
		start=end;
		
		logger.info("PCI-SPH integrate");
		CLEvent event = run_pcisph_integrate();
		end=System.currentTimeMillis();
		logger.info("PCI-SPH integrate end, took "+ (end-start) +"ms");
		start=end;
		
		// wait for the end of the run_pcisph_integrate on device
		event.waitFor();
		
		logger.info("SPH finish queue");
		// TODO: figure out if we need to actually call this
		_queue.finish();
		end=System.currentTimeMillis();
		logger.info("SPH finish queue end, took "+ (end-start) +"ms");
		start=end;
		
		endStep=System.currentTimeMillis();
		logger.info("SPH computation step done, took "+ (endStep-startStep) +"ms");
	}

	public void finishQueue()
	{
		_queue.finish();
	}

	public Float round(Float val, int roundingFactor)
	{
		return (float) Math.round(val * roundingFactor) / roundingFactor;
	}

	public int getParticleCountRoundedUp()
	{
		return (((_particleCount - 1) / 256) + 1) * 256;
	}

	@Override
	public StateSet solve(IRunConfiguration timeConfiguration)
	{
		// TODO: extend this to use time configuration to do multiple steps in one go
		long time = System.currentTimeMillis();
		logger.info("SPH solver start");
		StateSet stateSet = new StateSet(_model.getId());
		for(int i = 0; i < timeConfiguration.getTimeSteps(); i++)
		{
			long end=0;
			long start=System.currentTimeMillis();
			logger.info("SPH STEP START");
			step();
			updateStateSet(stateSet);
			end=System.currentTimeMillis();
			logger.info("SPH STEP END, took "+ (end-start) +"ms");
		}
		logger.info("SPH solver end, took: " + (System.currentTimeMillis() - time) + "ms");
		return stateSet;
	}

	private void updateStateSet(StateSet stateSet)
	{
		_positionPtr = _position.map(_queue, CLMem.MapFlags.Read);
		for(int i = 0, index=0; i < _particleCount; i++,index=index+4)
		{
			stateSet.addStateValue(new StateInstancePath(SPHModelInterpreterService.getPropertyPath(i, "pos", "x")), ValuesFactory.getFloatValue(_positionPtr.get(index)));
			stateSet.addStateValue(new StateInstancePath(SPHModelInterpreterService.getPropertyPath(i, "pos", "y")), ValuesFactory.getFloatValue(_positionPtr.get(index + 1)));
			stateSet.addStateValue(new StateInstancePath(SPHModelInterpreterService.getPropertyPath(i, "pos", "z")), ValuesFactory.getFloatValue(_positionPtr.get(index + 2)));
			stateSet.addStateValue(new StateInstancePath(SPHModelInterpreterService.getPropertyPath(i, "pos", "p")), ValuesFactory.getFloatValue(_positionPtr.get(index + 3)));
		}
		_position.unmap(_queue, _positionPtr);
	}

	@Override
	public void initialize(IModel model) throws GeppettoInitializationException
	{
		_model = (SPHModelX) model;
		setBuffersFromModel();
	}

	@Override
	public void dispose()
	{
		// TODO Erase all buffers close the context and buonanotte al secchio
		cleanContext();
	}
};