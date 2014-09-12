/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

package org.geppetto.solver.sph;

import static java.lang.System.out;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.StringTokenizer;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.bridj.Pointer;
import org.geppetto.core.common.GeppettoExecutionException;
import org.geppetto.core.common.GeppettoInitializationException;
import org.geppetto.core.data.model.AVariable;
import org.geppetto.core.data.model.ArrayVariable;
import org.geppetto.core.data.model.SimpleType;
import org.geppetto.core.data.model.SimpleType.Type;
import org.geppetto.core.data.model.SimpleVariable;
import org.geppetto.core.data.model.StructuredType;
import org.geppetto.core.data.model.VariableList;
import org.geppetto.core.model.IModel;
import org.geppetto.core.model.quantities.PhysicalQuantity;
import org.geppetto.core.model.runtime.ACompositeNode;
import org.geppetto.core.model.runtime.ANode;
import org.geppetto.core.model.runtime.AspectNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode;
import org.geppetto.core.model.runtime.CompositeNode;
import org.geppetto.core.model.runtime.ParticleNode;
import org.geppetto.core.model.runtime.VariableNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode.AspectTreeType;
import org.geppetto.core.model.values.FloatValue;
import org.geppetto.core.model.values.ValuesFactory;
import org.geppetto.core.simulation.IRunConfiguration;
import org.geppetto.core.solver.ISolver;
import org.geppetto.core.utilities.VariablePathSerializer;
import org.geppetto.core.visualisation.model.Point;
import org.geppetto.model.sph.Connection;
import org.geppetto.model.sph.Membrane;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.Vector3DX;
import org.springframework.stereotype.Service;
import org.springframework.web.servlet.view.velocity.VelocityConfig;

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
public class SPHSolverService implements ISolver {

	private static Log logger = LogFactory.getLog(SPHSolverService.class);

	private VariableList watchableVariables = new VariableList();
	private VariableList forceableVariables = new VariableList();

	List<String> watchListVarNames = new ArrayList<String>();
	boolean watching = false;

	private CLContext _context;
	public CLQueue _queue;
	private CLProgram _program;
	private CLDevice _device;
	private CLBuffer<Float> _acceleration;
	private CLBuffer<Integer> _gridCellIndex;
	private CLBuffer<Integer> _gridCellIndexFixedUp;
	private CLBuffer<Float> _neighborMap;
	private CLBuffer<Integer> _particleIndex;
	private CLBuffer<Integer> _particleIndexBack;
	private CLBuffer<Float> _position;
	private CLBuffer<Float> _pressure;
	private CLBuffer<Float> _rho;
	private CLBuffer<Float> _sortedPosition;
	private CLBuffer<Float> _sortedVelocity;
	private CLBuffer<Float> _velocity;
	private CLBuffer<Float> _elasticConnectionsData;
	private CLBuffer<Float> _activationSignal;
	private CLBuffer<Integer> _particleMembranesList; // potentially any particle can be connected with others via membrane(s)
													   // this buffer contains MAX_MEMBRANES_INCLUDING_SAME_PARTICLE integer data cells per particle
													   // each cell can contain -1 in case when no or no more membranes are associated with this particle,
													   // or the index of corresponding membrane in membraneData list othewize
	private CLBuffer<Integer> _membraneData;// elementary membrane is built on 3 adjacent particles (i,j,k) and should have a form of triangle
											 // highly recommended that i-j, j-k and k-i are already connected with springs to keep them close 
											 // to each other during whole lifetime of the simulation (user should control this by him(her)self)

	private Pointer<Float> _accelerationPtr;
	private Pointer<Integer> _gridCellIndexPtr;
	private Pointer<Integer> _gridCellIndexFixedUpPtr;
	private Pointer<Float> _neighborMapPtr;
	private Pointer<Integer> _particleIndexPtr;
	private Pointer<Integer> _particleIndexBackPtr;
	private Pointer<Float> _positionPtr;
	private Pointer<Float> _pressurePtr;
	private Pointer<Float> _rhoPtr;
	private Pointer<Float> _sortedPositionPtr;
	private Pointer<Float> _sortedVelocityPtr;
	private Pointer<Float> _velocityPtr;
	private Pointer<Float> _elasticConnectionsDataPtr;
	private Pointer<Float> _activationSignalPtr;
	private Pointer<Integer> _particleMembranesListPtr;
	private Pointer<Integer> _membraneDataPtr;

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
	
	//Membranes interaction kernels
	private CLKernel _clearMembraneBuffers;
	private CLKernel _computeInteractionWithMembranes;
	private CLKernel _computeInteractionWithMembranes_finalize;

	public float _xMax;
	public float _xMin;
	public float _yMax;
	public float _yMin;
	public float _zMax;
	public float _zMin;
	public int _elasticBundlesCount = 0;
	
	public float _mass;
	public float _timeStep;
	public float _simulationScale;
	public float _simulationScaleInv;
	public float _surfTensCoeff;
	public float _elasticityCoeff;
 	public float _viscosityCoeff;
	
	public int _gridCellsX;
	public int _gridCellsY;
	public int _gridCellsZ;
	public int _gridCellCount;
	public int _particleCount;
	public int _numOfLiquidP;
	public int _numOfElasticP;
	public int _numOfBoundaryP;
	public int _numOfMembranes;
	

	private SPHModelX _model;
	
	private int iterationNumber;
	private boolean _recordCheckPoints = false;
	private KernelsEnum checkKernel = null;
	private BuffersEnum checkBuffer = null;
	

	/*
	 * Checkpoints for the last computed step NOTE: stores all buffer values
	 * after each kernel execution for troubleshooting purposes
	 */
	private Map<KernelsEnum, PCISPHCheckPoint> _checkpointsMap = new LinkedHashMap<KernelsEnum, PCISPHCheckPoint>();

	public Map<KernelsEnum, PCISPHCheckPoint> getCheckpointsMap() {
		return _checkpointsMap;
	}
	public void setRecordCheckpoint(boolean value){
		this._recordCheckPoints = value;
	}
	/*
	 * A map of buffer sizes NOTE: these values are used in multiple places so
	 * storing them here reduces potential for error
	 */
	private Map<BuffersEnum, Integer> _buffersSizeMap = new LinkedHashMap<BuffersEnum, Integer>();

	public Map<BuffersEnum, Integer> getBuffersSizeMap() {
		return _buffersSizeMap;
	}

	public static Random RandomGenerator = new Random();

	public SPHSolverService(HardwareProfileEnum hardwareProfile)
			throws Exception {
		this.onceOffInit(hardwareProfile);
	}

	public SPHSolverService() throws Exception {
		this(HardwareProfileEnum.CPU);
	}

	public SPHSolverService(boolean recordCheckpoints) throws Exception {
		this();

		_recordCheckPoints = recordCheckpoints;
	}

	private void onceOffInit(HardwareProfileEnum hwProfile) throws IOException {
		// TODO: check if the selected profile is actually available
		DeviceFeature feature = (hwProfile == HardwareProfileEnum.CPU) ? DeviceFeature.CPU
				: DeviceFeature.CPU;

		_context = JavaCL.createBestContext(feature);

		out.println("created " + _context);
		// an array with available devices
		CLDevice[] devices = _context.getDevices();

		for (int i = 0; i < devices.length; i++) {
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
		String src = IOUtils.readText(SPHSolverService.class
				.getResourceAsStream("/resource/sphFluid.cl"));
		_program = _context.createProgram(src);

		// kernels
		_clearBuffers = _program.createKernel(KernelsEnum.CLEAR_BUFFERS
				.toString());
		_findNeighbors = _program.createKernel(KernelsEnum.FIND_NEIGHBORS
				.toString());
		_hashParticles = _program.createKernel(KernelsEnum.HASH_PARTICLES
				.toString());
		_indexx = _program.createKernel(KernelsEnum.INDEX.toString());
		_sortPostPass = _program.createKernel(KernelsEnum.SORT_POST_PASS
				.toString());

		// PCI-SPH specific
		_pcisph_computeForcesAndInitPressure = _program
				.createKernel(KernelsEnum.COMPUTE_FORCES_INIT_PRESSURE
						.toString());
		_pcisph_integrate = _program.createKernel(KernelsEnum.INTEGRATE
				.toString());
		_pcisph_predictPositions = _program
				.createKernel(KernelsEnum.PREDICT_POSITION.toString());
		_pcisph_predictDensity = _program
				.createKernel(KernelsEnum.PREDICT_DENSITY.toString());
		_pcisph_correctPressure = _program
				.createKernel(KernelsEnum.CORRECT_PRESSURE.toString());
		_pcisph_computePressureForceAcceleration = _program
				.createKernel(KernelsEnum.COMPUTE_PRESSURE_FORCE_ACCELERATION
						.toString());
		_pcisph_computeDensity = _program
				.createKernel(KernelsEnum.COMPUTE_DENSITY.toString());
		_pcisph_computeElasticForces = _program
				.createKernel(KernelsEnum.COMPUTE_ELASTIC_FORCES.toString());
		// Membranes interaction 
		_clearMembraneBuffers = _program.createKernel(KernelsEnum.CLEAR_MEMBRANE_BUFFERS.toString());
		_computeInteractionWithMembranes = _program.createKernel(KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES.toString());
		_computeInteractionWithMembranes_finalize = _program.createKernel( KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES_FINALIZE.toString() );
	}

	private void allocateBuffers() {
		// init buffer size map
		_buffersSizeMap.put(BuffersEnum.ACCELERATION, _particleCount * 4 * 3);
		_buffersSizeMap.put(BuffersEnum.GRID_CELL_INDEX, _gridCellCount + 1);
		_buffersSizeMap.put(BuffersEnum.GRID_CELL_INDEX_FIXED,
				_gridCellCount + 1);
		_buffersSizeMap.put(BuffersEnum.NEIGHBOR_MAP, _particleCount
				* SPHConstants.MAX_NEIGHBOR_COUNT * 2);
		_buffersSizeMap.put(BuffersEnum.PARTICLE_INDEX, _particleCount * 2);
		_buffersSizeMap.put(BuffersEnum.PARTICLE_INDEX_BACK, _particleCount);
		_buffersSizeMap.put(BuffersEnum.POSITION, _particleCount * 4 * ( 1 + 1 )/*1 extra, for membrane handling*/);
		_buffersSizeMap.put(BuffersEnum.PRESSURE, _particleCount * 1);
		_buffersSizeMap.put(BuffersEnum.RHO, _particleCount * 2);
		_buffersSizeMap
				.put(BuffersEnum.SORTED_POSITION, _particleCount * 4 * 2);
		_buffersSizeMap.put(BuffersEnum.SORTED_VELOCITY, _particleCount * 4);
		_buffersSizeMap.put(BuffersEnum.VELOCITY, _particleCount * 4 * ( 1 + 1 )/*1 extra, for membrane handling*/);
		_buffersSizeMap.put(BuffersEnum.ELASTIC_BUNDLES, _elasticBundlesCount);

		// allocate native device memory for all buffers
		_acceleration = _context.createFloatBuffer(CLMem.Usage.InputOutput,
				_buffersSizeMap.get(BuffersEnum.ACCELERATION));
		_gridCellIndex = _context.createIntBuffer(CLMem.Usage.InputOutput,
				_buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX));
		_gridCellIndexFixedUp = _context.createIntBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		_neighborMap = _context.createFloatBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.NEIGHBOR_MAP));
		_particleIndex = _context.createIntBuffer(CLMem.Usage.InputOutput,
				_buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX));
		_particleIndexBack = _context.createIntBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.PARTICLE_INDEX_BACK));
		_position = _context.createFloatBuffer(CLMem.Usage.InputOutput,
				_buffersSizeMap.get(BuffersEnum.POSITION));
		_pressure = _context.createFloatBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.PRESSURE));
		_rho = _context.createFloatBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.RHO));
		_sortedPosition = _context.createFloatBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.SORTED_POSITION));
		_sortedVelocity = _context.createFloatBuffer(
				_recordCheckPoints ? CLMem.Usage.InputOutput
						: CLMem.Usage.Input, _buffersSizeMap
						.get(BuffersEnum.SORTED_VELOCITY));
		_velocity = _context.createFloatBuffer(CLMem.Usage.InputOutput,
				_buffersSizeMap.get(BuffersEnum.VELOCITY));
	}

	private void setBuffersFromModel() {
		// set dimensions
		_xMax = _model.getXMax();
		_xMin = _model.getXMin();
		_yMax = _model.getYMax();
		_yMin = _model.getYMin();
		_zMax = _model.getZMax();
		_zMin = _model.getZMin();
		//Set simulation constant mass, timestep, simulationscle ... 
		_mass = (_model.getMass() != null) ? _model.getMass():SPHConstants.MASS;
		_timeStep = ( _model.getTimeStep() != null) ? _model.getTimeStep(): SPHConstants.TIME_STEP;
		_simulationScale = (float) ( 0.004f * Math.pow( _mass, 1.f/3.f ) / Math.pow( 0.00025f, 1.f/3.f ) );
		_simulationScaleInv = 1 / _simulationScale;
		_surfTensCoeff = ( _model.getSurfTensionCoeff() != null) ? _model.getSurfTensionCoeff(): SPHConstants.SURFACE_TENSION_COEFFICIENT;
		_elasticityCoeff = ( _model.getElasticitiCoeff() != null) ? _model.getElasticitiCoeff(): SPHConstants.ELASTICITY_COEFFICIENT;
		_elasticityCoeff /= _mass;
		_viscosityCoeff = ( _model.getViscosityCoeff() != null) ? _model.getViscosityCoeff(): SPHConstants.viscosity;
		_elasticBundlesCount = (_model.getElasticBundles() == null) ? 0
				: _model.getElasticBundles().intValue();
		SPHConstants.setBeta(_timeStep, _mass);
		SPHConstants.setDependingParammeters(_simulationScale, _mass);
		_surfTensCoeff *= (SPHConstants.W_POLY_6_COEFFICIENT * Math.pow((SPHConstants.H * SPHConstants.H * _simulationScale  * _simulationScale)/2.0,3.0)) * _simulationScale; 
		_particleCount = _model.getNumberOfParticles();
		_numOfElasticP = 0;
		_numOfLiquidP = 0;
		_numOfBoundaryP = 0;
		
		_gridCellsX = (int) ((_model.getXMax() - _model.getXMin()) / SPHConstants.H) + 1;
		_gridCellsY = (int) ((_model.getYMax() - _model.getYMin()) / SPHConstants.H) + 1;
		_gridCellsZ = (int) ((_model.getZMax() - _model.getZMin()) / SPHConstants.H) + 1;

		// set grid dimensions
		_gridCellCount = _gridCellsX * _gridCellsY * _gridCellsZ;
		
		//OUT all main constant
		System.out.println("/******MAIN CONSTAN FOR CONFIGURATION************/");
		System.out.println("particle mass: " + _mass);
		System.out.println("time step: " + _timeStep);
		System.out.println("simulation scale: " + _simulationScale);
		System.out.println("surface tension coefficient: " + _surfTensCoeff);
		System.out.println("elasticity coefficient: " + _elasticityCoeff);
		System.out.println("viscosity coefficient: " + _viscosityCoeff);
		System.out.println("beta: " + SPHConstants.BETA);
		System.out.println("delta: " + SPHConstants.DELTA);
		System.out.println("W_POLY_6_COEFFICIENT: " + SPHConstants.W_POLY_6_COEFFICIENT);
		System.out.println("GRAD_W_SPIKY_COEFFICIENT: " + SPHConstants.GRAD_W_SPIKY_COEFFICIENT);
		System.out.println("DEL_2_W_VISCOSITY_COEFFICIENT: " + SPHConstants.DEL_2_W_VISCOSITY_COEFFICIENT);
		System.out.println("MASS_MULT_WPOLY6COEFFICIENT: " + SPHConstants.MASS_MULT_WPOLY6COEFFICIENT);
		System.out.println("MASS_MULT_GRADWSPIKYCOEFFICIENT: " + SPHConstants.MASS_MULT_GRADWSPIKYCOEFFICIENT);
		System.out.println("MASS_MULT_DIVGRADWVISCOSITYCOEFFICIENT: " + SPHConstants.MASS_MULT_DIVGRADWVISCOSITYCOEFFICIENT);
		System.out.println("/*******************END**************************/");
		// allocate buffers - requires global dimensions of the grid
		this.allocateBuffers();

		int index = 0;

		for (int i = 0; i < _particleCount; i++) {
			if (i != 0) {
				index = index + 4;
			}

			Vector3DX positionVector = (Vector3DX) _model.getParticles().get(i)
					.getPositionVector();
			Vector3DX velocityVector = (Vector3DX) _model.getParticles().get(i)
					.getVelocityVector();

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
			if (positionVector.getP() == SPHConstants.BOUNDARY_TYPE) {
				_numOfBoundaryP++;
			} else if (positionVector.getP() == SPHConstants.ELASTIC_TYPE) {
				_numOfElasticP++;
			} else if (positionVector.getP() == SPHConstants.LIQUID_TYPE) {
				_numOfLiquidP++;
			}
		}

		// populate elastic connection buffers if we have any
		if (_numOfElasticP > 0 && _model.getConnections().size() > 0) {
			// init elastic connections buffers
			// TODO: move this back with the other buffers init stuff
			_buffersSizeMap.put(BuffersEnum.ELASTIC_CONNECTIONS, _numOfElasticP
					* SPHConstants.MAX_NEIGHBOR_COUNT * 4);
			_elasticConnectionsData = _context.createFloatBuffer(
					CLMem.Usage.InputOutput,
					_buffersSizeMap.get(BuffersEnum.ELASTIC_CONNECTIONS));
			_elasticConnectionsDataPtr = _elasticConnectionsData.map(_queue,
					CLMem.MapFlags.Write);

			int connIndex = 0;
			for (Connection conn : _model.getConnections()) {
				_elasticConnectionsDataPtr.set(connIndex, conn.getP1());
				_elasticConnectionsDataPtr.set(connIndex + 1,
						conn.getDistance());
				_elasticConnectionsDataPtr.set(connIndex + 2,
						conn.getMysteryValue());
				_elasticConnectionsDataPtr.set(connIndex + 3, 0f); // padding
				connIndex += 4;
			}
			//init Actication signal buffer
			_buffersSizeMap.put(BuffersEnum.ACTIVATION_SIGNAL, SPHConstants.MUSCLE_COUNT);
			_activationSignal = _context.createFloatBuffer(
					CLMem.Usage.InputOutput,
					_buffersSizeMap.get(BuffersEnum.ACTIVATION_SIGNAL));
			_activationSignalPtr = _activationSignal.map(_queue,
					CLMem.MapFlags.Write);
			for(int i=0;i<SPHConstants.MUSCLE_COUNT;i++){
				_activationSignalPtr.set(0.f);
			}
			// we copied the stuff down to the device and we won't touch it
			// again so we can unmap
			_elasticConnectionsData.unmap(_queue, _elasticConnectionsDataPtr);

			// allocate activation signal buffers
			if (_buffersSizeMap.get(BuffersEnum.ELASTIC_BUNDLES) > 0) {
				_activationSignal = _context.createFloatBuffer(
						CLMem.Usage.Input,
						_buffersSizeMap.get(BuffersEnum.ELASTIC_BUNDLES));
			} else {
				// allocate a buffer with 1 single element
				// NOTE: this is a HACK to avoid exceptions in case of having
				// elastic particles but no contractible bundles
				_activationSignal = _context.createFloatBuffer(
						CLMem.Usage.Input, 1);
			}
			if( _model.getMembranes().size() > 0 ) {
				//init membranes data buffer
				_numOfMembranes = _model.getMembranes().size();
				_buffersSizeMap.put(BuffersEnum.MEMBRANES_DATA, _numOfMembranes * 3);
				_membraneData = _context.createIntBuffer(
						CLMem.Usage.InputOutput,
						_buffersSizeMap.get(BuffersEnum.MEMBRANES_DATA));
				_membraneDataPtr = _membraneData.map(_queue,
						CLMem.MapFlags.Write);
				int memIndex = 0;
				for (Membrane m: _model.getMembranes()){
					_membraneDataPtr.set(memIndex, m.getParticleI());
					_membraneDataPtr.set(memIndex + 1,
							m.getParticleJ());
					_membraneDataPtr.set(memIndex + 2,
							m.getParticleK());
					memIndex += 3;
				}
				_membraneData.unmap(_queue, _membraneDataPtr);
	
				//init membranes particle buffer it contains indexes of membranes which contains particle
				_buffersSizeMap.put(BuffersEnum.MEMBRANES_PARTICLE_INDEX_LIST, _numOfElasticP
						* SPHConstants.MAX_MEMBRANES_INCLUDING_SAME_PARTICLE);
				_particleMembranesList = _context.createIntBuffer(
						CLMem.Usage.InputOutput,
						_buffersSizeMap.get(BuffersEnum.MEMBRANES_PARTICLE_INDEX_LIST));
				_particleMembranesListPtr = _particleMembranesList.map(_queue,
						CLMem.MapFlags.Write);
				int _index = 0;
				for (Integer mIndex : _model.getParticleMembranesList()) {
					_particleMembranesListPtr.set(_index,mIndex);
					++_index;
				}
				_particleMembranesList.unmap(_queue, _particleMembranesListPtr);
			}
		}

		// check that counts are fine
		if (_particleCount != (_numOfBoundaryP + _numOfElasticP + _numOfLiquidP)) {
			throw new IllegalArgumentException(
					"SPHSolverService:setModels - particle counts do not add up");
		}
	}

	public void cleanContext() {
		_context.release();
	}
	//TODO depricated method will be removed after determenistic test
	private int runClearBuffers() {
		
		_clearBuffers.setArg(0, _neighborMap);
		_clearBuffers.setArg(1, _particleCount);
		_clearBuffers.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	private CLEvent runHashParticles() {
		// Stage HashParticles
		_hashParticles.setArg(0, _position);
		_hashParticles.setArg(1, _gridCellsX);
		_hashParticles.setArg(2, _gridCellsY);
		_hashParticles.setArg(3, _gridCellsZ);
		_hashParticles.setArg(4, SPHConstants.HASH_GRID_CELL_SIZE_INV);
		_hashParticles.setArg(5, _xMin);
		_hashParticles.setArg(6, _yMin);
		_hashParticles.setArg(7, _zMin);
		_hashParticles.setArg(8, _particleIndex);
		_hashParticles.setArg(9, _particleCount);
		CLEvent event = _hashParticles.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return event;
	}
	private CLEvent runIndexx() {
		// Stage Indexx
		_indexx.setArg(0, _particleIndex);
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_indexx.setArg(1, _gridCellCount);
		_indexx.setArg(2, _gridCellIndex);
		_indexx.setArg(3, _particleCount);
		int gridCellCountRoundedUp = (((_gridCellCount - 1) / 256) + 1) * 256;
		CLEvent event = _indexx.enqueueNDRange(_queue,
				new int[] { gridCellCountRoundedUp });

		return event;
	}
	private int runIndexPostPass() {
		// get values out of buffer
		_gridCellIndexPtr = _gridCellIndex.map(_queue, CLMem.MapFlags.Read);
		int[] gridNextNonEmptyCellBuffer = _gridCellIndexPtr.getInts();
		_gridCellIndex.unmap(_queue, _gridCellIndexPtr);

		int recentNonEmptyCell = _gridCellCount;
		for (int i = _gridCellCount; i >= 0; i--) {
			if (gridNextNonEmptyCellBuffer[i] == SPHConstants.NO_CELL_ID) {
				gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell;
			} else {
				recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
			}
		}

		// put results back
		_gridCellIndexFixedUpPtr = _gridCellIndexFixedUp.map(_queue,
				CLMem.MapFlags.Write);
		_gridCellIndexFixedUpPtr.setInts(gridNextNonEmptyCellBuffer);
		_gridCellIndexFixedUp.unmap(_queue, _gridCellIndexFixedUpPtr);

		return 0;
	}
	
	private int runSort() {
		// this version work with qsort
		int index = 0;
		List<int[]> particleIndex = new ArrayList<int[]>();

		// get values out of buffer
		_particleIndexPtr = _particleIndex
				.map(_queue, CLMem.MapFlags.ReadWrite);
		int[] particleInd = _particleIndexPtr.getInts();

		for (int i = 0; i < _particleCount * 2; i += 2) {
			int[] element = { particleInd[i], particleInd[i + 1] };
			particleIndex.add(element);
		}
		Collections.sort(particleIndex, new MyCompare());
		for (int i = 0; i < particleIndex.size(); i++) {
			for (int j = 0; j < 2; j++) {
				particleInd[index] = particleIndex.get(i)[j];
				index++;
			}
		}

		// put results back
		_particleIndexPtr.setInts(particleInd);
		_particleIndex.unmap(_queue, _particleIndexPtr);

		return 0;
	}
	
	private int runSortPostPass() {
		// Stage SortPostPass
		_sortPostPass.setArg(0, _particleIndex);
		_sortPostPass.setArg(1, _particleIndexBack);
		_sortPostPass.setArg(2, _position);
		_sortPostPass.setArg(3, _velocity);
		_sortPostPass.setArg(4, _sortedPosition);
		_sortPostPass.setArg(5, _sortedVelocity);
		_sortPostPass.setArg(6, _particleCount);
		_sortPostPass.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });
		return 0;
	}
	
	private int runFindNeighbors() {
		_findNeighbors.setArg(0, _gridCellIndexFixedUp);
		_findNeighbors.setArg(1, _sortedPosition);
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_findNeighbors.setArg(2, _gridCellCount);
		_findNeighbors.setArg(3, _gridCellsX);
		_findNeighbors.setArg(4, _gridCellsY);
		_findNeighbors.setArg(5, _gridCellsZ);
		_findNeighbors.setArg(6, SPHConstants.H);
		_findNeighbors.setArg(7, SPHConstants.HASH_GRID_CELL_SIZE);
		_findNeighbors.setArg(8, SPHConstants.HASH_GRID_CELL_SIZE_INV);
		_findNeighbors.setArg(9, _simulationScale);
		_findNeighbors.setArg(10, _xMin);
		_findNeighbors.setArg(11, _yMin);
		_findNeighbors.setArg(12, _zMin);
		_findNeighbors.setArg(13, _neighborMap);
		_findNeighbors.setArg(14, _particleCount);
		_findNeighbors.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });
		return 0;
	}
	
	private int run_pcisph_computeDensity() {
		// Stage ComputeDensityPressure
		_pcisph_computeDensity.setArg(0, _neighborMap);
		_pcisph_computeDensity.setArg(1, SPHConstants.MASS_MULT_WPOLY6COEFFICIENT);
		_pcisph_computeDensity.setArg(2, SPHConstants._hScaled2);
		_pcisph_computeDensity.setArg(3, _rho);
		_pcisph_computeDensity.setArg(4, _particleIndexBack);
		_pcisph_computeDensity.setArg(5, _particleCount);
		_pcisph_computeDensity.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computeForcesAndInitPressure() {
		_pcisph_computeForcesAndInitPressure.setArg(0, _neighborMap);
		_pcisph_computeForcesAndInitPressure.setArg(1, _rho);
		_pcisph_computeForcesAndInitPressure.setArg(2, _pressure);
		_pcisph_computeForcesAndInitPressure.setArg(3, _sortedPosition);
		_pcisph_computeForcesAndInitPressure.setArg(4, _sortedVelocity);
		_pcisph_computeForcesAndInitPressure.setArg(5, _acceleration);
		_pcisph_computeForcesAndInitPressure.setArg(6, _particleIndexBack);
		_pcisph_computeForcesAndInitPressure.setArg(7, _surfTensCoeff);
		_pcisph_computeForcesAndInitPressure.setArg(8,
				SPHConstants.MASS_MULT_DIVGRADWVISCOSITYCOEFFICIENT);
		_pcisph_computeForcesAndInitPressure.setArg(9, SPHConstants._hScaled);
		_pcisph_computeForcesAndInitPressure.setArg(10, _viscosityCoeff);
		_pcisph_computeForcesAndInitPressure.setArg(11, SPHConstants.GRAVITY_X);
		_pcisph_computeForcesAndInitPressure.setArg(12, SPHConstants.GRAVITY_Y);
		_pcisph_computeForcesAndInitPressure.setArg(13, SPHConstants.GRAVITY_Z);
		_pcisph_computeForcesAndInitPressure.setArg(14, _position);
		_pcisph_computeForcesAndInitPressure.setArg(15, _particleIndex);
		_pcisph_computeForcesAndInitPressure.setArg(16, _particleCount);
		_pcisph_computeForcesAndInitPressure.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computeElasticForces() {
		_pcisph_computeElasticForces.setArg(0, _neighborMap);
		_pcisph_computeElasticForces.setArg(1, _sortedPosition);
		_pcisph_computeElasticForces.setArg(2, _sortedVelocity);
		_pcisph_computeElasticForces.setArg(3, _acceleration);
		_pcisph_computeElasticForces.setArg(4, _particleIndexBack);
		_pcisph_computeElasticForces.setArg(5, _velocity);
		_pcisph_computeElasticForces.setArg(6, SPHConstants.H);
		_pcisph_computeElasticForces.setArg(7, _mass);
		_pcisph_computeElasticForces.setArg(8, _simulationScale);
		_pcisph_computeElasticForces.setArg(9, _numOfElasticP);
		_pcisph_computeElasticForces.setArg(10, _elasticConnectionsData);
		_pcisph_computeElasticForces.setArg(11, _particleCount);
		_pcisph_computeElasticForces.setArg(12, SPHConstants.MUSCLE_COUNT);
		_pcisph_computeElasticForces.setArg(13, _activationSignal);
		_pcisph_computeElasticForces.setArg(14, _position);
		_pcisph_computeElasticForces.setArg(15, _elasticityCoeff);
		
		int numOfElasticPRoundedUp = (((_numOfElasticP - 1) / 256) + 1) * 256;

		_pcisph_computeElasticForces.enqueueNDRange(_queue,
				new int[] { numOfElasticPRoundedUp });

		return 0;
	}

	private int run_pcisph_predictPositions() {
		_pcisph_predictPositions.setArg(0, _acceleration);
		_pcisph_predictPositions.setArg(1, _sortedPosition);
		_pcisph_predictPositions.setArg(2, _sortedVelocity);
		_pcisph_predictPositions.setArg(3, _particleIndex);
		_pcisph_predictPositions.setArg(4, _particleIndexBack);
		_pcisph_predictPositions.setArg(5, SPHConstants.GRAVITY_X);
		_pcisph_predictPositions.setArg(6, SPHConstants.GRAVITY_Y);
		_pcisph_predictPositions.setArg(7, SPHConstants.GRAVITY_Z);
		_pcisph_predictPositions.setArg(8, this._simulationScaleInv);
		_pcisph_predictPositions.setArg(9, this._timeStep);
		_pcisph_predictPositions.setArg(10, _position);
		_pcisph_predictPositions.setArg(11, _velocity);
		_pcisph_predictPositions.setArg(12, SPHConstants.R0);
		_pcisph_predictPositions.setArg(13, _neighborMap);
		_pcisph_predictPositions.setArg(14, _particleCount);
		_pcisph_predictPositions.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_predictDensity() {
		// Stage predict density
		_pcisph_predictDensity.setArg(0, _neighborMap);
		_pcisph_predictDensity.setArg(1, _particleIndexBack);
		_pcisph_predictDensity.setArg(2, SPHConstants.MASS_MULT_WPOLY6COEFFICIENT);
		_pcisph_predictDensity.setArg(3, SPHConstants.H);
		_pcisph_predictDensity.setArg(4, SPHConstants.RHO0);
		_pcisph_predictDensity.setArg(5, this._simulationScale);
		_pcisph_predictDensity.setArg(6, SPHConstants.STIFFNESS);
		_pcisph_predictDensity.setArg(7, _sortedPosition);
		_pcisph_predictDensity.setArg(8, _pressure);
		_pcisph_predictDensity.setArg(9, _rho);
		_pcisph_predictDensity.setArg(10, this._particleCount);
		_pcisph_predictDensity.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_correctPressure() {
		// Stage correct pressure
		_pcisph_correctPressure.setArg(0, _particleIndexBack);
		_pcisph_correctPressure.setArg(1, SPHConstants.RHO0);
		_pcisph_correctPressure.setArg(2, _pressure);
		_pcisph_correctPressure.setArg(3, _rho);
		_pcisph_correctPressure.setArg(4, SPHConstants.DELTA);
		_pcisph_correctPressure.setArg(5, _particleCount);
		_pcisph_correctPressure.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computePressureForceAcceleration() {
		// Stage ComputeAcceleration
		_pcisph_computePressureForceAcceleration.setArg(0, _neighborMap);
		_pcisph_computePressureForceAcceleration.setArg(1, _pressure);
		_pcisph_computePressureForceAcceleration.setArg(2, _rho);
		_pcisph_computePressureForceAcceleration.setArg(3, _sortedPosition);
		_pcisph_computePressureForceAcceleration.setArg(4, _sortedVelocity);
		_pcisph_computePressureForceAcceleration.setArg(5, _particleIndexBack);
		_pcisph_computePressureForceAcceleration.setArg(6,
				SPHConstants.DELTA);
		_pcisph_computePressureForceAcceleration.setArg(7,
				SPHConstants.MASS_MULT_GRADWSPIKYCOEFFICIENT);
		_pcisph_computePressureForceAcceleration.setArg(8, SPHConstants.H);
		_pcisph_computePressureForceAcceleration.setArg(9, _simulationScale);
		_pcisph_computePressureForceAcceleration.setArg(10, this._viscosityCoeff);
		_pcisph_computePressureForceAcceleration.setArg(11, _acceleration);
		_pcisph_computePressureForceAcceleration.setArg(12, SPHConstants.RHO0);
		_pcisph_computePressureForceAcceleration.setArg(13, _position);
		_pcisph_computePressureForceAcceleration.setArg(14, _particleIndex);
		_pcisph_computePressureForceAcceleration.setArg(15, _particleCount);
		_pcisph_computePressureForceAcceleration.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private CLEvent run_pcisph_integrate() {
		// Stage Integrate
		_pcisph_integrate.setArg(0, _acceleration);
		_pcisph_integrate.setArg(1, _sortedPosition);
		_pcisph_integrate.setArg(2, _sortedVelocity);
		_pcisph_integrate.setArg(3, _particleIndex);
		_pcisph_integrate.setArg(4, _particleIndexBack);
		_pcisph_integrate.setArg(5, SPHConstants.GRAVITY_X);
		_pcisph_integrate.setArg(6, SPHConstants.GRAVITY_Y);
		_pcisph_integrate.setArg(7, SPHConstants.GRAVITY_Z);
		_pcisph_integrate.setArg(8, this._simulationScaleInv);
		_pcisph_integrate.setArg(9, this._timeStep);
		_pcisph_integrate.setArg(10, _xMin);
		_pcisph_integrate.setArg(11, _xMax);
		_pcisph_integrate.setArg(12, _yMin);
		_pcisph_integrate.setArg(13, _yMax);
		_pcisph_integrate.setArg(14, _zMin);
		_pcisph_integrate.setArg(15, _zMax);
		_pcisph_integrate.setArg(16, SPHConstants.DAMPING);
		_pcisph_integrate.setArg(17, _position);
		_pcisph_integrate.setArg(18, _velocity);
		_pcisph_integrate.setArg(19, _rho);
		_pcisph_integrate.setArg(20, SPHConstants.R0);
		_pcisph_integrate.setArg(21, _neighborMap);
		_pcisph_integrate.setArg(22, _particleCount);
		_pcisph_integrate.setArg(23, iterationNumber);
		CLEvent event = _pcisph_integrate.enqueueNDRange(_queue,
				new int[] { getParticleCountRoundedUp() });

		return event;
	}
	private int run_clearMembraneBuffers(){
		// Stage Clear membrane Buffeers
		_clearMembraneBuffers.setArg(0, _position);
		_clearMembraneBuffers.setArg(1, _velocity);
		_clearMembraneBuffers.setArg(2, _sortedPosition);
		_clearMembraneBuffers.setArg(3, _particleCount);
		_clearMembraneBuffers.enqueueNDRange( _queue, new int[]{ getParticleCountRoundedUp() } );
		return 0;
	}
	private int run_computeInteractionWithMembranes(){
		_computeInteractionWithMembranes.setArg( 0, _position );
		_computeInteractionWithMembranes.setArg( 1, _velocity );
		_computeInteractionWithMembranes.setArg( 2, _sortedPosition );
		_computeInteractionWithMembranes.setArg( 3, _particleIndex );
		_computeInteractionWithMembranes.setArg( 4, _particleIndexBack );
		_computeInteractionWithMembranes.setArg( 5, _neighborMap );
		_computeInteractionWithMembranes.setArg( 6, _particleMembranesList );
		_computeInteractionWithMembranes.setArg( 7, _membraneData );
		_computeInteractionWithMembranes.setArg( 8, _particleCount );
		_computeInteractionWithMembranes.setArg( 9, _numOfElasticP );
		_computeInteractionWithMembranes.setArg( 10, SPHConstants.R0 );
		_computeInteractionWithMembranes.enqueueNDRange(_queue, new int[]{ getParticleCountRoundedUp() } );
		return 0;
	}
	private CLEvent run_computeInteractionWithMembranes_finalize(){
		_computeInteractionWithMembranes_finalize.setArg( 0, _position );
		_computeInteractionWithMembranes_finalize.setArg( 1, _velocity );
		_computeInteractionWithMembranes_finalize.setArg( 2, _particleIndex );
		_computeInteractionWithMembranes_finalize.setArg( 3, _particleIndexBack );
		_computeInteractionWithMembranes_finalize.setArg( 4, _particleCount );
		CLEvent event = _computeInteractionWithMembranes_finalize.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return event;
	}
	
	class MyCompare implements Comparator<int[]> {
		public int compare(int[] o1, int[] o2) {
			if (o1[0] < o2[0])
				return -1;
			if (o1[0] > o2[0])
				return +1;
			return 0;
		}
	}

	private void step() {
		long endStep = 0;
		long startStep = System.currentTimeMillis();
		long end = 0;
		long start = System.currentTimeMillis();
		//DEPRICATED FUNCTION
		/*logger.info("SPH clear buffer");
		runClearBuffers();
		if (_recordCheckPoints) {
			recordCheckpoints(KernelsEnum.CLEAR_BUFFERS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH clear buffer end, took " + (end - start) + "ms");
		start = end;*/

		logger.info("SPH hash particles");
		CLEvent hashParticles = runHashParticles();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.HASH_PARTICLES)
				recordCheckpoints(KernelsEnum.HASH_PARTICLES);
		}
		end = System.currentTimeMillis();
		logger.info("SPH hash particles end, took " + (end - start) + "ms");
		start = end;

		// host needs to wait as the next operation requires values from buffers
		hashParticles.waitFor();

		logger.info("SPH sort");
		runSort();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.SORT)
				recordCheckpoints(KernelsEnum.SORT);
		}
		end = System.currentTimeMillis();
		logger.info("SPH sort end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH sort post pass");
		runSortPostPass();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.SORT_POST_PASS)
				recordCheckpoints(KernelsEnum.SORT_POST_PASS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH sort post pass end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH index");
		CLEvent runIndexx = runIndexx();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.INDEX)
				recordCheckpoints(KernelsEnum.INDEX);
		}
		end = System.currentTimeMillis();
		logger.info("SPH index end, took " + (end - start) + "ms");
		start = end;

		// host needs to wait as the next operation requires values from buffers
		runIndexx.waitFor();

		logger.info("SPH index post pass");
		runIndexPostPass();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.INDEX_POST_PASS)
				recordCheckpoints(KernelsEnum.INDEX_POST_PASS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH index post pass end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH find neighbors");
		runFindNeighbors();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.FIND_NEIGHBORS)
				recordCheckpoints(KernelsEnum.FIND_NEIGHBORS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH find neighbors end, took " + (end - start) + "ms");
		start = end;

		// PCISPH stuff starts here
		logger.info("PCI-SPH compute density");
		run_pcisph_computeDensity();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.COMPUTE_DENSITY)
				recordCheckpoints(KernelsEnum.COMPUTE_DENSITY);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH compute density end, took " + (end - start) + "ms");
		start = end;

		logger.info("PCI-SPH compute forces and init pressure");
		run_pcisph_computeForcesAndInitPressure();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.COMPUTE_FORCES_INIT_PRESSURE)
				recordCheckpoints(KernelsEnum.COMPUTE_FORCES_INIT_PRESSURE);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH compute forces and init pressure end, took "
				+ (end - start) + "ms");
		start = end;

		// Do elastic stuff only if we have elastic particles
		if (_numOfElasticP > 0) {
			logger.info("PCI-SPH compute elastic forces");
			run_pcisph_computeElasticForces();
			if (_recordCheckPoints) {
				if(checkKernel == null || checkKernel == KernelsEnum.COMPUTE_ELASTIC_FORCES)
					recordCheckpoints(KernelsEnum.COMPUTE_ELASTIC_FORCES);
			}
			end = System.currentTimeMillis();
			logger.info("PCI-SPH compute elastic forces end, took "
					+ (end - start) + "ms");
			start = end;
		}

		logger.info("PCI-SPH predict/correct loop");
		// LOOP: 3 times or until "error" becomes less than 2%
		int iter = 0;
		int maxIterations = 3;
		do {
			run_pcisph_predictPositions();
			run_pcisph_predictDensity();
			run_pcisph_correctPressure();
			run_pcisph_computePressureForceAcceleration();

			iter++;
		} while ((iter < maxIterations));
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.PREDICTIVE_LOOP)
				recordCheckpoints(KernelsEnum.PREDICTIVE_LOOP);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH predict/correct loop end, took " + (end - start)
				+ "ms");
		start = end;

		logger.info("PCI-SPH integrate");
		CLEvent event = run_pcisph_integrate();
		if (_recordCheckPoints) {
			if(checkKernel == null || checkKernel == KernelsEnum.INTEGRATE)
				recordCheckpoints(KernelsEnum.INTEGRATE);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH integrate end, took " + (end - start) + "ms");
		start = end;
		// wait for the end of the run_pcisph_integrate on device
		event.waitFor();
		if(_numOfMembranes!=0){
		   logger.info("PCI-SPH membrane interaction calculating");
		   run_clearMembraneBuffers();
		   if (_recordCheckPoints) {
			   if(checkKernel == null || checkKernel == KernelsEnum.CLEAR_MEMBRANE_BUFFERS)
				   recordCheckpoints(KernelsEnum.CLEAR_MEMBRANE_BUFFERS);
		   }
		   run_computeInteractionWithMembranes();
		   if (_recordCheckPoints) {
			   if(checkKernel == null || checkKernel == KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES)
				   recordCheckpoints(KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES);
		   }
		   // compute change of coordinates due to interactions with membranes
		   event = run_computeInteractionWithMembranes_finalize();
		   if (_recordCheckPoints) {
			   if(checkKernel == null || checkKernel == KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES_FINALIZE)
				   recordCheckpoints(KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES_FINALIZE);
		   }
		   // wait for the end of the run_pcisph_integrate on device
		   event.waitFor();
		   end = System.currentTimeMillis();
		   logger.info("PCI-SPH membrane intraction calculating end, took " + (end - start)
					+ "ms");
		   start = end;
		}

		logger.info("SPH finish queue");
		// TODO: figure out if we need to actually call this
		_queue.finish();
		end = System.currentTimeMillis();
		logger.info("SPH finish queue end, took " + (end - start) + "ms");
		start = end;

		endStep = System.currentTimeMillis();
		logger.info("SPH computation step done, took " + (endStep - startStep)
				+ "ms");
	}

	public void finishQueue() {
		_queue.finish();
	}

	private int getParticleCountRoundedUp() {
		return (((_particleCount - 1) / 256) + 1) * 256;
	}

	@Override
	public void solve(IRunConfiguration timeConfiguration, AspectNode aspect) throws GeppettoExecutionException {
		// TODO: extend this to use time configuration to do multiple steps in one go
		long time = System.currentTimeMillis();
		logger.info("SPH solver start");

		for (iterationNumber = 0; iterationNumber < timeConfiguration.getTimeSteps(); iterationNumber++) {
			// TODO: setActivationSignal

			long end = 0;
			long start = System.currentTimeMillis();
			logger.info("SPH STEP START");
			step();
			
			updateStateTree(aspect);

			end = System.currentTimeMillis();
			logger.info("Step count " + iterationNumber);
			logger.info("SPH STEP END, took " + (end - start) + "ms");
			
			
		}

		logger.info("SPH solver end, took: " + (System.currentTimeMillis() - time) + "ms");
	}

	private void updateStateTree(AspectNode aspect) {
		AspectSubTreeNode visualTree = (AspectSubTreeNode) aspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		AspectSubTreeNode simulationTree = (AspectSubTreeNode) aspect.getSubTree(AspectTreeType.WATCH_TREE);

		_positionPtr = _position.map(_queue, CLMem.MapFlags.Read);

		// ASSUMPTION: The solver will never create new states after the first
		// time step
		// we can call it principle of conservation of the states; if there is a
		// good
		// reason to revoke this assumption we need to add code that at every
		// cycle checks
		// if some new states exist to eventually add them to the stateTree
		UpdateSPHVisualizationTreeVisitor updateSPHStateTreeVisitor = new UpdateSPHVisualizationTreeVisitor(_positionPtr);
		visualTree.apply(updateSPHStateTreeVisitor);

		if (watching) {
			updateSimulationTree(simulationTree);
		}
		
		_position.unmap(_queue, _positionPtr);

		visualTree.setModified(true);
		simulationTree.setModified(true);
	}

	/**
	 * Updates nodes for simulation tree, if tree is empty call method to populate it.
	 * 
	 * @param simulationTree
	 */
	private void updateSimulationTree(AspectSubTreeNode simulationTree) {
		// map watchable buffers that are not already mapped
		// NOTE: position is mapped for scene generation - improving performance by not mapping it again
		_velocityPtr = _velocity.map(_queue, CLMem.MapFlags.Read);

		if (simulationTree.getChildren().isEmpty()) {	
			populateSimulationTree(simulationTree);			
		} else {
			// watch tree not empty populate new values
			UpdateSPHSimulationTreeVisitor visitor = new UpdateSPHSimulationTreeVisitor(_positionPtr);
			simulationTree.apply(visitor);
		}

		// unmap watchable buffers
		_velocity.unmap(_queue, _positionPtr);
	}

	private boolean containsNode(ACompositeNode node, String name){
		List<ANode> children = node.getChildren();
		
		boolean addNewNode = true;
		for(ANode child : children){
			if(child.getName().equals(name)){
				addNewNode = false;
				return addNewNode;
			}
			if(child instanceof ACompositeNode){
				if(((ACompositeNode)child).getChildren() != null){
					addNewNode = containsNode((ACompositeNode) child, name);
				}
			}

		}
		
		return addNewNode;
	}
	
	private ACompositeNode getNode(ACompositeNode node, String name){
		ACompositeNode newNode = null;
		
		List<ANode> children = node.getChildren();
		
		boolean addNewNode = true;
		for(ANode child : children){
			if(child.getName().equals(name)){
				newNode = (ACompositeNode) child;
				return newNode;
			}
			if(child instanceof ACompositeNode){
				if(((ACompositeNode)child).getChildren() != null){
					newNode = getNode((ACompositeNode) child, name);
				}
			}

		}
		
		return newNode;
	}
	
	@Override
	public void initialize(IModel model) {
		_model = (SPHModelX) model;
		
		setBuffersFromModel();

		setWatchableVariables();
		setForceableVariables();

	}

	@Override
	public void dispose() {
		// close the context and "buonanotte al secchio" (good night to the
		// bucket)
		cleanContext();
	}

	private void setActivationSignal(float[] activation) {
		// put results back
		_activationSignalPtr = _activationSignal.map(_queue,
				CLMem.MapFlags.Write);
		_activationSignalPtr.setFloats(activation);
		_activationSignal.unmap(_queue, _activationSignalPtr);
	}

	private void recordCheckpoints(KernelsEnum kernelCheckpoint) {
		PCISPHCheckPoint check = new PCISPHCheckPoint();

		// read buffers into lists and populate checkpoint object
		if(checkBuffer == BuffersEnum.ACCELERATION || checkBuffer == null)
			check.acceleration = this.<Float> getBufferValues(_accelerationPtr,
					_acceleration,
					this._buffersSizeMap.get(BuffersEnum.ACCELERATION));
		if(checkBuffer == BuffersEnum.GRID_CELL_INDEX || checkBuffer == null)
			check.gridCellIndex = this.<Integer> getBufferValues(_gridCellIndexPtr,
					_gridCellIndex,
					this._buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX));
		if(checkBuffer ==BuffersEnum.GRID_CELL_INDEX_FIXED || checkBuffer == null)
			check.gridCellIndexFixedUp = this.<Integer> getBufferValues(
					_gridCellIndexFixedUpPtr, _gridCellIndexFixedUp,
					this._buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		if(checkBuffer ==BuffersEnum.NEIGHBOR_MAP || checkBuffer == null)
			check.neighborMap = this.<Float> getBufferValues(_neighborMapPtr,
					_neighborMap,
					this._buffersSizeMap.get(BuffersEnum.NEIGHBOR_MAP));
		if(checkBuffer ==BuffersEnum.PARTICLE_INDEX || checkBuffer == null)
			check.particleIndex = this.<Integer> getBufferValues(_particleIndexPtr,
					_particleIndex,
					this._buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX));
		if(checkBuffer ==BuffersEnum.PARTICLE_INDEX_BACK || checkBuffer == null)
			check.particleIndexBack = this.<Integer> getBufferValues(
					_particleIndexBackPtr, _particleIndexBack,
					this._buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX_BACK));
		if(checkBuffer ==BuffersEnum.POSITION || checkBuffer == null)
			check.position = this.<Float> getBufferValues(_positionPtr, _position,
					this._buffersSizeMap.get(BuffersEnum.POSITION));
		if(checkBuffer ==BuffersEnum.PRESSURE || checkBuffer == null)
			check.pressure = this.<Float> getBufferValues(_pressurePtr, _pressure,
					this._buffersSizeMap.get(BuffersEnum.PRESSURE));
		if(checkBuffer ==BuffersEnum.RHO || checkBuffer == null)
			check.rho = this.<Float> getBufferValues(_rhoPtr, _rho,
					this._buffersSizeMap.get(BuffersEnum.RHO));
		if(checkBuffer ==BuffersEnum.SORTED_POSITION || checkBuffer == null)
			check.sortedPosition = this.<Float> getBufferValues(_sortedPositionPtr,
					_sortedPosition,
					this._buffersSizeMap.get(BuffersEnum.SORTED_POSITION));
		if(checkBuffer ==BuffersEnum.SORTED_VELOCITY || checkBuffer == null)
			check.sortedVelocity = this.<Float> getBufferValues(_sortedVelocityPtr,
					_sortedVelocity,
					this._buffersSizeMap.get(BuffersEnum.SORTED_VELOCITY));
		if(checkBuffer ==BuffersEnum.VELOCITY || checkBuffer == null)
			check.velocity = this.<Float> getBufferValues(_velocityPtr, _velocity,
					this._buffersSizeMap.get(BuffersEnum.VELOCITY));
		if (_numOfElasticP > 0) {
			if(checkBuffer ==BuffersEnum.ELASTIC_CONNECTIONS || checkBuffer == null)
				check.elasticConnections = this.<Float> getBufferValues(
						_elasticConnectionsDataPtr, _elasticConnectionsData,
						this._buffersSizeMap.get(BuffersEnum.ELASTIC_CONNECTIONS));
		}
		if(_numOfMembranes > 0){
			if(checkBuffer ==BuffersEnum.MEMBRANES_DATA || checkBuffer == null)
				check.membranes = this.<Integer> getBufferValues(_membraneDataPtr, _membraneData, 
						this._buffersSizeMap.get(BuffersEnum.MEMBRANES_DATA));
			if(checkBuffer ==BuffersEnum.MEMBRANES_PARTICLE_INDEX_LIST || checkBuffer == null)
				check.membranesParticleIndexList = this.<Integer> getBufferValues(_particleMembranesListPtr, _particleMembranesList, 
						this._buffersSizeMap.get(BuffersEnum.MEMBRANES_PARTICLE_INDEX_LIST));
		}
		_checkpointsMap.put(kernelCheckpoint, check);
	}

	/*
	 * A method to retrieve buffer values into simple lists
	 */
	private <T> List<T> getBufferValues(Pointer<T> pointer, CLBuffer<T> buffer,
			int size) {
		List<T> list = new ArrayList<T>();

		pointer = buffer.map(_queue, CLMem.MapFlags.Read);

		for (int i = 0; i < size; i++) {
			list.add(pointer.get(i));
		}

		buffer.unmap(_queue, pointer);

		return list;
	}

	@Override
	public VariableList getForceableVariables() {
		return forceableVariables;
	}

	@Override
	public VariableList getWatchableVariables() {
		return watchableVariables;
	}

	/**
	 * Populates state variables that can be watched
	 * 
	 * */
	private void setWatchableVariables() {

		SimpleType floatType = new SimpleType();
		floatType.setType(Type.FLOAT);

		// structure type vector
		StructuredType vector = new StructuredType();
		List<AVariable> vectorVars = new ArrayList<AVariable>();
		SimpleVariable x = new SimpleVariable();
		SimpleVariable y = new SimpleVariable();
		SimpleVariable z = new SimpleVariable();
		x.setName("x");
		x.setType(floatType);
		y.setName("y");
		y.setType(floatType);
		z.setName("z");
		z.setType(floatType);
		vectorVars.addAll(Arrays.asList(x, y, z));
		vector.setVariables(vectorVars);

		// structure type particle
		StructuredType particle = new StructuredType();
		List<AVariable> particleVars = new ArrayList<AVariable>();
		SimpleVariable position = new SimpleVariable();
		SimpleVariable velocity = new SimpleVariable();
		position.setName("position");
		position.setType(vector);
		velocity.setName("velocity");
		velocity.setType(vector);
		particleVars.addAll(Arrays.asList(position, velocity));
		particle.setVariables(particleVars);

		List<AVariable> vars = new ArrayList<AVariable>();

		// array of particles
		ArrayVariable particles = new ArrayVariable();
		particles.setName("particle");
		particles.setType(particle);
		particles.setSize(_particleCount);

		vars.add(particles);

		this.watchableVariables.setVariables(vars);
	}

	/**
	 * Populates state variables that can be watched
	 * 
	 * */
	private void setForceableVariables() {
		List<AVariable> vars = new ArrayList<AVariable>();

		// a float type
		SimpleType floatType = new SimpleType();
		floatType.setType(Type.FLOAT);

		// activation signals - array of floats
		ArrayVariable activationSignals = new ArrayVariable();
		activationSignals.setName("activation");
		activationSignals.setType(floatType);
		activationSignals.setSize(_buffersSizeMap
				.get(BuffersEnum.ELASTIC_BUNDLES));

		vars.add(activationSignals);

		this.forceableVariables.setVariables(vars);
	}

	@Override
	public void addWatchVariables(List<String> variableNames) {
		watchListVarNames.addAll(variableNames);
	}

	@Override
	public void startWatch() {
		watching = true;
	}

	@Override
	public void stopWatch() {
		watching = false;
	}

	@Override
	public void clearWatchVariables() {
		watchListVarNames.clear();
	}
	/**
	 * @return the checkKernel
	 */
	public KernelsEnum getCheckKernel() {
		return checkKernel;
	}
	/**
	 * @param checkKernel the checkKernel to set
	 */
	public void setCheckKernel(KernelsEnum checkKernel) {
		this.checkKernel = checkKernel;
	}
	/**
	 * @return the checkBuffer
	 */
	public BuffersEnum getCheckBuffer() {
		return checkBuffer;
	}
	/**
	 * @param checkBuffer the checkBuffer to set
	 */
	public void setCheckBuffer(BuffersEnum checkBuffer) {
		this.checkBuffer = checkBuffer;
	}

	@Override
	public void populateVisualTree(IModel model, AspectSubTreeNode visualTree) throws GeppettoInitializationException
	{
		_positionPtr = _position.map(_queue, CLMem.MapFlags.Read);

		CompositeNode _liquidModel = new CompositeNode("LIQUID_" + model.getId());
		CompositeNode _boundaryModel = new CompositeNode("BOUNDARY_" + model.getId());
		CompositeNode _elasticModel = new CompositeNode("ELASTIC_" + model.getId());
		
		// the state tree is empty, let's create it
		for (int i = 0, index = 0; i < _particleCount; i++, index = index + 4) {
			String particleId = SPHModelInterpreterService.getParticleId(i);
			FloatValue xV = ValuesFactory.getFloatValue(_positionPtr.get(index));
			FloatValue yV = ValuesFactory.getFloatValue(_positionPtr.get(index + 1));
			FloatValue zV = ValuesFactory.getFloatValue(_positionPtr.get(index + 2));
			FloatValue pV = ValuesFactory.getFloatValue(_positionPtr.get(index + 3));

			if (pV.getAsFloat() != SPHConstants.BOUNDARY_TYPE) {
				// don't need to create a state for the boundary particles,
				// they don't move.
				ParticleNode particle = new ParticleNode(particleId);
				Point pos=new Point();
				pos.setX(xV.getAsDouble());
				pos.setY(yV.getAsDouble());
				pos.setZ(zV.getAsDouble());
				particle.setPosition(pos);
				particle.setParticleKind(pV.getAsFloat());
				particle.setId(particleId);
				
				if(pV.getAsFloat() == (SPHConstants.LIQUID_TYPE))
				{
					_liquidModel.addChild(particle);
				}
				else if(pV.getAsFloat() == (SPHConstants.ELASTIC_TYPE))
				{
					_elasticModel.addChild(particle);
				}
				else if(pV.getAsFloat() == (SPHConstants.BOUNDARY_TYPE))
				{
					_boundaryModel.addChild(particle);
				}
			}
		} 
		
		visualTree.addChild(_liquidModel);
		visualTree.addChild(_elasticModel);
		visualTree.addChild(_boundaryModel);		
		
		_position.unmap(_queue, _positionPtr);
	}

	@Override
	public void populateSimulationTree(AspectSubTreeNode watchTree)
	{
		// map watchable buffers that are not already mapped
		// NOTE: position is mapped for scene generation - improving performance by not mapping it again
		_velocityPtr = _velocity.map(_queue, CLMem.MapFlags.Read);

		// check which watchable variables are being watched
		for (AVariable var : getWatchableVariables().getVariables()) {
			for (String varName : watchListVarNames) {
				// get watchable variables path
				List<String> watchableVarsPaths = new ArrayList<String>();
				VariablePathSerializer.GetFullVariablePath(var, "", watchableVarsPaths);

				varName = varName.replace(watchTree.getInstancePath()+".", "");
				// remove array bracket arguments from variable paths
				String varNameNoBrackets = varName;
				String particleID = null;
				if(varName.indexOf("[")!=-1)
				{
					varNameNoBrackets = varName.substring(0,varName.indexOf("["))+varName.substring(varName.indexOf("]")+1,varName.length());
					particleID = varName.substring(varName.indexOf("[")+1, varName.indexOf("]"));
				}

				// loop through paths and look for matching paths
				for (String s : watchableVarsPaths) {
					if (s.equals(varNameNoBrackets)) {
						// we have a match

						Integer ID = null;
						if (particleID != null) {
							// check that paticleID is valid
							ID = Integer.parseInt(particleID);
							if (!(ID < _particleCount)) {
								throw new IllegalArgumentException("SPHSolverService:updateStateTreeForWatch - particle index is out of boundaries");
							}
						}

						// tokenize variable path in watch list via dot
						// separator (handle array brackets)
						StringTokenizer tokenizer = new StringTokenizer(s,".");
						ACompositeNode node = watchTree;
						while (tokenizer.hasMoreElements()) {
							// loop through tokens and build tree
							String current = tokenizer.nextToken();
							boolean found = false;
							for (ANode child : node.getChildren()) {
								if (child.getName().equals(current)) {
									if (child instanceof ACompositeNode) {
										node = (ACompositeNode) child;
									}
									found = true;
									break;
								}
							}
							if (found) {
								continue;
							} else {
								if (tokenizer.hasMoreElements()) {
									// not a leaf, create a composite statenode
									String nodeName = current;
									if(current.equals("particle"))
									{
										nodeName = current + "[" + particleID + "]";
									}

									CompositeNode newNode = new CompositeNode(nodeName);
									newNode.setId(nodeName);
									
									boolean addNewNode = containsNode(node, newNode.getName());

									if(addNewNode){
										node.addChild(newNode);
										node = newNode;
									}
									else{
										node = getNode(node, newNode.getName());
									}
								} else {
									// it's a leaf node
									VariableNode newNode = new VariableNode(current);
									newNode.setId(current);
									
									FloatValue val = null;

									// get value
									switch (current) {
										case "x":
											val = ValuesFactory.getFloatValue(_positionPtr.get(ID));
											break;
										case "y":
											val = ValuesFactory.getFloatValue(_positionPtr.get(ID + 1));
											break;
										case "z":
											val = ValuesFactory.getFloatValue(_positionPtr.get(ID + 2));
											break;
									}

									PhysicalQuantity q = new PhysicalQuantity();
									q.setValue(val);
									newNode.addPhysicalQuantity(q);

									node.addChild(newNode);
								}
							}
						}
					}
				}
			}
		}
		
		watchTree.setModified(true);
		// unmap watchable buffers
		_velocity.unmap(_queue, _positionPtr);
	}
};
