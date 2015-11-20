/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011 - 2015 OpenWorm.
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
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.bridj.Pointer;
import org.geppetto.core.common.GeppettoExecutionException;
import org.geppetto.core.data.model.IAspectConfiguration;
import org.geppetto.core.model.IModel;
import org.geppetto.core.model.runtime.ANode;
import org.geppetto.core.model.runtime.AspectSubTreeNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode.AspectTreeType;
import org.geppetto.core.model.runtime.EntityNode;
import org.geppetto.core.model.typesystem.AspectNode;
import org.geppetto.core.model.typesystem.values.ACompositeValue;
import org.geppetto.core.services.AService;
import org.geppetto.core.solver.ISolver;
import org.geppetto.model.sph.Connection;
import org.geppetto.model.sph.common.SPHConstants;
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
public class SPHSolverService extends AService implements ISolver
{

	private static Log logger = LogFactory.getLog(SPHSolverService.class);

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
	public int _elasticBundlesCount = 0;

	public int _gridCellsX;
	public int _gridCellsY;
	public int _gridCellsZ;
	public int _gridCellCount;
	public int _particleCount;
	public int _numOfLiquidP;
	public int _numOfElasticP;
	public int _numOfBoundaryP;

	private SPHModelX _model;

	private boolean _recordCheckPoints = false;

	/*
	 * Checkpoints for the last computed step NOTE: stores all buffer values after each kernel execution for troubleshooting purposes
	 */
	private Map<KernelsEnum, PCISPHCheckPoint> _checkpointsMap = new LinkedHashMap<KernelsEnum, PCISPHCheckPoint>();

	public Map<KernelsEnum, PCISPHCheckPoint> getCheckpointsMap()
	{
		return _checkpointsMap;
	}

	/*
	 * A map of buffer sizes NOTE: these values are used in multiple places so storing them here reduces potential for error
	 */
	private Map<BuffersEnum, Integer> _buffersSizeMap = new LinkedHashMap<BuffersEnum, Integer>();

	public Map<BuffersEnum, Integer> getBuffersSizeMap()
	{
		return _buffersSizeMap;
	}

	public static Random RandomGenerator = new Random();

	public SPHSolverService(HardwareProfileEnum hardwareProfile) throws Exception
	{
		this.onceOffInit(hardwareProfile);
	}

	public SPHSolverService() throws Exception
	{
		this(HardwareProfileEnum.GPU);
	}

	public SPHSolverService(boolean recordCheckpoints) throws Exception
	{
		this();

		_recordCheckPoints = recordCheckpoints;
	}

	private void onceOffInit(HardwareProfileEnum hwProfile) throws IOException
	{
		// TODO: check if the selected profile is actually available
		DeviceFeature feature = (hwProfile == HardwareProfileEnum.CPU) ? DeviceFeature.CPU : DeviceFeature.GPU;

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
		_clearBuffers = _program.createKernel(KernelsEnum.CLEAR_BUFFERS.toString());
		_findNeighbors = _program.createKernel(KernelsEnum.FIND_NEIGHBORS.toString());
		_hashParticles = _program.createKernel(KernelsEnum.HASH_PARTICLES.toString());
		_indexx = _program.createKernel(KernelsEnum.INDEX.toString());
		_sortPostPass = _program.createKernel(KernelsEnum.SORT_POST_PASS.toString());

		// PCI-SPH specific
		_pcisph_computeForcesAndInitPressure = _program.createKernel(KernelsEnum.COMPUTE_FORCES_INIT_PRESSURE.toString());
		_pcisph_integrate = _program.createKernel(KernelsEnum.INTEGRATE.toString());
		_pcisph_predictPositions = _program.createKernel(KernelsEnum.PREDICT_POSITION.toString());
		_pcisph_predictDensity = _program.createKernel(KernelsEnum.PREDICT_DENSITY.toString());
		_pcisph_correctPressure = _program.createKernel(KernelsEnum.CORRECT_PRESSURE.toString());
		_pcisph_computePressureForceAcceleration = _program.createKernel(KernelsEnum.COMPUTE_PRESSURE_FORCE_ACCELERATION.toString());
		_pcisph_computeDensity = _program.createKernel(KernelsEnum.COMPUTE_DENSITY.toString());
		_pcisph_computeElasticForces = _program.createKernel(KernelsEnum.COMPUTE_ELASTIC_FORCES.toString());
	}

	private void allocateBuffers()
	{
		// init buffer size map
		_buffersSizeMap.put(BuffersEnum.ACCELERATION, _particleCount * 4 * 2);
		_buffersSizeMap.put(BuffersEnum.GRID_CELL_INDEX, _gridCellCount + 1);
		_buffersSizeMap.put(BuffersEnum.GRID_CELL_INDEX_FIXED, _gridCellCount + 1);
		_buffersSizeMap.put(BuffersEnum.NEIGHBOR_MAP, _particleCount * SPHConstants.NEIGHBOR_COUNT * 2);
		_buffersSizeMap.put(BuffersEnum.PARTICLE_INDEX, _particleCount * 2);
		_buffersSizeMap.put(BuffersEnum.PARTICLE_INDEX_BACK, _particleCount);
		_buffersSizeMap.put(BuffersEnum.POSITION, _particleCount * 4);
		_buffersSizeMap.put(BuffersEnum.PRESSURE, _particleCount * 4);
		_buffersSizeMap.put(BuffersEnum.RHO, _particleCount * 2);
		_buffersSizeMap.put(BuffersEnum.SORTED_POSITION, _particleCount * 4 * 2);
		_buffersSizeMap.put(BuffersEnum.SORTED_VELOCITY, _particleCount * 4);
		_buffersSizeMap.put(BuffersEnum.VELOCITY, _particleCount * 4);
		_buffersSizeMap.put(BuffersEnum.ELASTIC_BUNDLES, _elasticBundlesCount);

		// allocate native device memory for all buffers
		_acceleration = _context.createFloatBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.ACCELERATION));
		_gridCellIndex = _context.createIntBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX));
		_gridCellIndexFixedUp = _context.createIntBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		_neighborMap = _context.createFloatBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.NEIGHBOR_MAP));
		_particleIndex = _context.createIntBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX));
		_particleIndexBack = _context.createIntBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX_BACK));
		_position = _context.createFloatBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.POSITION));
		_pressure = _context.createFloatBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.PRESSURE));
		_rho = _context.createFloatBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.RHO));
		_sortedPosition = _context.createFloatBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.SORTED_POSITION));
		_sortedVelocity = _context.createFloatBuffer(_recordCheckPoints ? CLMem.Usage.InputOutput : CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.SORTED_VELOCITY));
		_velocity = _context.createFloatBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.VELOCITY));
	}

	private void setBuffersFromModel()
	{
		// set dimensions
		_xMax = _model.getXMax();
		_xMin = _model.getXMin();
		_yMax = _model.getYMax();
		_yMin = _model.getYMin();
		_zMax = _model.getZMax();
		_zMin = _model.getZMin();
		_elasticBundlesCount = (_model.getElasticBundles() == null) ? 0 : _model.getElasticBundles().intValue();

		_particleCount = _model.getNumberOfParticles();
		_numOfElasticP = 0;
		_numOfLiquidP = 0;
		_numOfBoundaryP = 0;

		_gridCellsX = (int) ((_model.getXMax() - _model.getXMin()) / SPHConstants.H) + 1;
		_gridCellsY = (int) ((_model.getYMax() - _model.getYMin()) / SPHConstants.H) + 1;
		_gridCellsZ = (int) ((_model.getZMax() - _model.getZMin()) / SPHConstants.H) + 1;

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
			_buffersSizeMap.put(BuffersEnum.ELASTIC_CONNECTIONS, _numOfElasticP * SPHConstants.NEIGHBOR_COUNT * 4);
			_elasticConnectionsData = _context.createFloatBuffer(CLMem.Usage.InputOutput, _buffersSizeMap.get(BuffersEnum.ELASTIC_CONNECTIONS));
			_elasticConnectionsDataPtr = _elasticConnectionsData.map(_queue, CLMem.MapFlags.Write);

			int connIndex = 0;
			for(Connection conn : _model.getConnections())
			{
				_elasticConnectionsDataPtr.set(connIndex, conn.getP1());
				_elasticConnectionsDataPtr.set(connIndex + 1, conn.getDistance());
				_elasticConnectionsDataPtr.set(connIndex + 2, conn.getMysteryValue());
				_elasticConnectionsDataPtr.set(connIndex + 3, 0f); // padding
				connIndex += 4;
			}

			// we copied the stuff down to the device and we won't touch it
			// again so we can unmap
			_elasticConnectionsData.unmap(_queue, _elasticConnectionsDataPtr);

			// allocate activation signal buffers
			if(_buffersSizeMap.get(BuffersEnum.ELASTIC_BUNDLES) > 0)
			{
				_activationSignal = _context.createFloatBuffer(CLMem.Usage.Input, _buffersSizeMap.get(BuffersEnum.ELASTIC_BUNDLES));
			}
			else
			{
				// allocate a buffer with 1 single element
				// NOTE: this is a HACK to avoid exceptions in case of having
				// elastic particles but no contractible bundles
				_activationSignal = _context.createFloatBuffer(CLMem.Usage.Input, 1);
			}
		}

		// check that counts are fine
		if(_particleCount != (_numOfBoundaryP + _numOfElasticP + _numOfLiquidP))
		{
			throw new IllegalArgumentException("SPHSolverService:setModels - particle counts do not add up");
		}
	}

	public void cleanContext()
	{
		_context.release();
	}

	private int runClearBuffers()
	{
		_clearBuffers.setArg(0, _neighborMap);
		_clearBuffers.setArg(1, _particleCount);
		_clearBuffers.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	private int runFindNeighbors()
	{
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
		_findNeighbors.setArg(9, SPHConstants.SIMULATION_SCALE);
		_findNeighbors.setArg(10, _xMin);
		_findNeighbors.setArg(11, _yMin);
		_findNeighbors.setArg(12, _zMin);
		_findNeighbors.setArg(13, _neighborMap);
		_findNeighbors.setArg(14, _particleCount);
		_findNeighbors.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });
		return 0;
	}

	private CLEvent runHashParticles()
	{
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
		CLEvent event = _hashParticles.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return event;
	}

	private int runIndexPostPass()
	{
		// get values out of buffer
		_gridCellIndexPtr = _gridCellIndex.map(_queue, CLMem.MapFlags.Read);
		int[] gridNextNonEmptyCellBuffer = _gridCellIndexPtr.getInts();
		_gridCellIndex.unmap(_queue, _gridCellIndexPtr);

		int recentNonEmptyCell = _gridCellCount;
		for(int i = _gridCellCount; i >= 0; i--)
		{
			if(gridNextNonEmptyCellBuffer[i] == SPHConstants.NO_CELL_ID)
			{
				gridNextNonEmptyCellBuffer[i] = recentNonEmptyCell;
			}
			else
			{
				recentNonEmptyCell = gridNextNonEmptyCellBuffer[i];
			}
		}

		// put results back
		_gridCellIndexFixedUpPtr = _gridCellIndexFixedUp.map(_queue, CLMem.MapFlags.Write);
		_gridCellIndexFixedUpPtr.setInts(gridNextNonEmptyCellBuffer);
		_gridCellIndexFixedUp.unmap(_queue, _gridCellIndexFixedUpPtr);

		return 0;
	}

	private CLEvent runIndexx()
	{
		// Stage Indexx
		_indexx.setArg(0, _particleIndex);
		_gridCellCount = ((_gridCellsX) * (_gridCellsY)) * (_gridCellsZ);
		_indexx.setArg(1, _gridCellCount);
		_indexx.setArg(2, _gridCellIndex);
		_indexx.setArg(3, _particleCount);
		int gridCellCountRoundedUp = (((_gridCellCount - 1) / 256) + 1) * 256;
		CLEvent event = _indexx.enqueueNDRange(_queue, new int[] { gridCellCountRoundedUp });

		return event;
	}

	private int runSortPostPass()
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

	private int run_pcisph_computeDensity()
	{
		// Stage ComputeDensityPressure
		_pcisph_computeDensity.setArg(0, _neighborMap);
		_pcisph_computeDensity.setArg(1, SPHConstants.W_POLY_6_COEFFICIENT);
		_pcisph_computeDensity.setArg(2, SPHConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_computeDensity.setArg(3, SPHConstants.H);
		_pcisph_computeDensity.setArg(4, SPHConstants.MASS);
		_pcisph_computeDensity.setArg(5, SPHConstants.RHO0);
		_pcisph_computeDensity.setArg(6, SPHConstants.SIMULATION_SCALE);
		_pcisph_computeDensity.setArg(7, SPHConstants.STIFFNESS);
		_pcisph_computeDensity.setArg(8, _sortedPosition);
		_pcisph_computeDensity.setArg(9, _pressure);
		_pcisph_computeDensity.setArg(10, _rho);
		_pcisph_computeDensity.setArg(11, _particleIndexBack);
		_pcisph_computeDensity.setArg(12, SPHConstants.DELTA); // calculated
																// from
																// constants
		_pcisph_computeDensity.setArg(13, _particleCount);
		_pcisph_computeDensity.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computeForcesAndInitPressure()
	{
		_pcisph_computeForcesAndInitPressure.setArg(0, _neighborMap);
		_pcisph_computeForcesAndInitPressure.setArg(1, _rho);
		_pcisph_computeForcesAndInitPressure.setArg(2, _pressure);
		_pcisph_computeForcesAndInitPressure.setArg(3, _sortedPosition);
		_pcisph_computeForcesAndInitPressure.setArg(4, _sortedVelocity);
		_pcisph_computeForcesAndInitPressure.setArg(5, _acceleration);
		_pcisph_computeForcesAndInitPressure.setArg(6, _particleIndexBack);
		_pcisph_computeForcesAndInitPressure.setArg(7, SPHConstants.W_POLY_6_COEFFICIENT);
		_pcisph_computeForcesAndInitPressure.setArg(8, SPHConstants.DEL_2_W_VISCOSITY_COEFFICIENT);
		_pcisph_computeForcesAndInitPressure.setArg(9, SPHConstants.H);
		_pcisph_computeForcesAndInitPressure.setArg(10, SPHConstants.MASS);
		_pcisph_computeForcesAndInitPressure.setArg(11, SPHConstants.MU);
		_pcisph_computeForcesAndInitPressure.setArg(12, SPHConstants.SIMULATION_SCALE);
		_pcisph_computeForcesAndInitPressure.setArg(13, SPHConstants.GRAVITY_X);
		_pcisph_computeForcesAndInitPressure.setArg(14, SPHConstants.GRAVITY_Y);
		_pcisph_computeForcesAndInitPressure.setArg(15, SPHConstants.GRAVITY_Z);
		_pcisph_computeForcesAndInitPressure.setArg(16, _position);
		_pcisph_computeForcesAndInitPressure.setArg(17, _particleIndex);
		_pcisph_computeForcesAndInitPressure.setArg(18, _particleCount);
		_pcisph_computeForcesAndInitPressure.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computeElasticForces()
	{
		_pcisph_computeElasticForces.setArg(0, _neighborMap);
		_pcisph_computeElasticForces.setArg(1, _sortedPosition);
		_pcisph_computeElasticForces.setArg(2, _sortedVelocity);
		_pcisph_computeElasticForces.setArg(3, _acceleration);
		_pcisph_computeElasticForces.setArg(4, _particleIndexBack);
		_pcisph_computeElasticForces.setArg(5, _velocity);
		_pcisph_computeElasticForces.setArg(6, SPHConstants.H);
		_pcisph_computeElasticForces.setArg(7, SPHConstants.MASS);
		_pcisph_computeElasticForces.setArg(8, SPHConstants.SIMULATION_SCALE);
		_pcisph_computeElasticForces.setArg(9, _numOfElasticP);
		_pcisph_computeElasticForces.setArg(10, _elasticConnectionsData);
		_pcisph_computeElasticForces.setArg(11, 0);
		_pcisph_computeElasticForces.setArg(12, _activationSignal);
		_pcisph_computeElasticForces.setArg(13, _elasticBundlesCount);
		_pcisph_computeElasticForces.setArg(14, _particleCount);

		int numOfElasticPRoundedUp = (((_numOfElasticP - 1) / 256) + 1) * 256;

		_pcisph_computeElasticForces.enqueueNDRange(_queue, new int[] { numOfElasticPRoundedUp });

		return 0;
	}

	private int run_pcisph_predictPositions()
	{
		_pcisph_predictPositions.setArg(0, _acceleration);
		_pcisph_predictPositions.setArg(1, _sortedPosition);
		_pcisph_predictPositions.setArg(2, _sortedVelocity);
		_pcisph_predictPositions.setArg(3, _particleIndex);
		_pcisph_predictPositions.setArg(4, _particleIndexBack);
		_pcisph_predictPositions.setArg(5, SPHConstants.GRAVITY_X);
		_pcisph_predictPositions.setArg(6, SPHConstants.GRAVITY_Y);
		_pcisph_predictPositions.setArg(7, SPHConstants.GRAVITY_Z);
		_pcisph_predictPositions.setArg(8, SPHConstants.SIMULATION_SCALE_INV);
		_pcisph_predictPositions.setArg(9, SPHConstants.TIME_STEP);
		_pcisph_predictPositions.setArg(10, _xMin);
		_pcisph_predictPositions.setArg(11, _xMax);
		_pcisph_predictPositions.setArg(12, _yMin);
		_pcisph_predictPositions.setArg(13, _yMax);
		_pcisph_predictPositions.setArg(14, _zMin);
		_pcisph_predictPositions.setArg(15, _zMax);
		_pcisph_predictPositions.setArg(16, SPHConstants.DAMPING);
		_pcisph_predictPositions.setArg(17, _position);
		_pcisph_predictPositions.setArg(18, _velocity);
		_pcisph_predictPositions.setArg(19, SPHConstants.R0);
		_pcisph_predictPositions.setArg(20, _neighborMap);
		_pcisph_predictPositions.setArg(21, _particleCount);
		_pcisph_predictPositions.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_predictDensity()
	{
		// Stage predict density
		_pcisph_predictDensity.setArg(0, _neighborMap);
		_pcisph_predictDensity.setArg(1, _particleIndexBack);
		_pcisph_predictDensity.setArg(2, SPHConstants.W_POLY_6_COEFFICIENT);
		_pcisph_predictDensity.setArg(3, SPHConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_predictDensity.setArg(4, SPHConstants.H);
		_pcisph_predictDensity.setArg(5, SPHConstants.MASS);
		_pcisph_predictDensity.setArg(6, SPHConstants.RHO0);
		_pcisph_predictDensity.setArg(7, SPHConstants.SIMULATION_SCALE);
		_pcisph_predictDensity.setArg(8, SPHConstants.STIFFNESS);
		_pcisph_predictDensity.setArg(9, _sortedPosition);
		_pcisph_predictDensity.setArg(10, _pressure);
		_pcisph_predictDensity.setArg(11, _rho);
		_pcisph_predictDensity.setArg(12, SPHConstants.DELTA);
		_pcisph_predictDensity.setArg(13, _particleCount);
		_pcisph_predictDensity.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_correctPressure()
	{
		// Stage correct pressure
		_pcisph_correctPressure.setArg(0, _neighborMap);
		_pcisph_correctPressure.setArg(1, _particleIndexBack);
		_pcisph_correctPressure.setArg(2, SPHConstants.W_POLY_6_COEFFICIENT);
		_pcisph_correctPressure.setArg(3, SPHConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_correctPressure.setArg(4, SPHConstants.H);
		_pcisph_correctPressure.setArg(5, SPHConstants.MASS);
		_pcisph_correctPressure.setArg(6, SPHConstants.RHO0);
		_pcisph_correctPressure.setArg(7, SPHConstants.SIMULATION_SCALE);
		_pcisph_correctPressure.setArg(8, SPHConstants.STIFFNESS);
		_pcisph_correctPressure.setArg(9, _sortedPosition);
		_pcisph_correctPressure.setArg(10, _pressure);
		_pcisph_correctPressure.setArg(11, _rho);
		_pcisph_correctPressure.setArg(12, SPHConstants.DELTA);
		_pcisph_correctPressure.setArg(13, _position);
		_pcisph_correctPressure.setArg(14, _particleIndex);
		_pcisph_correctPressure.setArg(15, _particleCount);
		_pcisph_correctPressure.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private int run_pcisph_computePressureForceAcceleration()
	{
		// Stage ComputeAcceleration
		_pcisph_computePressureForceAcceleration.setArg(0, _neighborMap);
		_pcisph_computePressureForceAcceleration.setArg(1, _pressure);
		_pcisph_computePressureForceAcceleration.setArg(2, _rho);
		_pcisph_computePressureForceAcceleration.setArg(3, _sortedPosition);
		_pcisph_computePressureForceAcceleration.setArg(4, _sortedVelocity);
		_pcisph_computePressureForceAcceleration.setArg(5, _particleIndexBack);
		_pcisph_computePressureForceAcceleration.setArg(6, SPHConstants.CFLLimit);
		_pcisph_computePressureForceAcceleration.setArg(7, SPHConstants.DEL_2_W_VISCOSITY_COEFFICIENT);
		_pcisph_computePressureForceAcceleration.setArg(8, SPHConstants.GRAD_W_SPIKY_COEFFICIENT);
		_pcisph_computePressureForceAcceleration.setArg(9, SPHConstants.H);
		_pcisph_computePressureForceAcceleration.setArg(10, SPHConstants.MASS);
		_pcisph_computePressureForceAcceleration.setArg(11, SPHConstants.MU);
		_pcisph_computePressureForceAcceleration.setArg(12, SPHConstants.SIMULATION_SCALE);
		_pcisph_computePressureForceAcceleration.setArg(13, _acceleration);
		_pcisph_computePressureForceAcceleration.setArg(14, SPHConstants.RHO0);
		_pcisph_computePressureForceAcceleration.setArg(15, _position);
		_pcisph_computePressureForceAcceleration.setArg(16, _particleIndex);
		_pcisph_computePressureForceAcceleration.setArg(17, _particleCount);
		_pcisph_computePressureForceAcceleration.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return 0;
	}

	private CLEvent run_pcisph_integrate()
	{
		// Stage Integrate
		_pcisph_integrate.setArg(0, _acceleration);
		_pcisph_integrate.setArg(1, _sortedPosition);
		_pcisph_integrate.setArg(2, _sortedVelocity);
		_pcisph_integrate.setArg(3, _particleIndex);
		_pcisph_integrate.setArg(4, _particleIndexBack);
		_pcisph_integrate.setArg(5, SPHConstants.GRAVITY_X);
		_pcisph_integrate.setArg(6, SPHConstants.GRAVITY_Y);
		_pcisph_integrate.setArg(7, SPHConstants.GRAVITY_Z);
		_pcisph_integrate.setArg(8, SPHConstants.SIMULATION_SCALE_INV);
		_pcisph_integrate.setArg(9, SPHConstants.TIME_STEP);
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
		CLEvent event = _pcisph_integrate.enqueueNDRange(_queue, new int[] { getParticleCountRoundedUp() });

		return event;
	}

	private int runSort()
	{
		// this version work with qsort
		int index = 0;
		List<int[]> particleIndex = new ArrayList<int[]>();

		// get values out of buffer
		_particleIndexPtr = _particleIndex.map(_queue, CLMem.MapFlags.ReadWrite);
		int[] particleInd = _particleIndexPtr.getInts();

		for(int i = 0; i < _particleCount * 2; i += 2)
		{
			int[] element = { particleInd[i], particleInd[i + 1] };
			particleIndex.add(element);
		}
		Collections.sort(particleIndex, new MyCompare());
		for(int i = 0; i < particleIndex.size(); i++)
		{
			for(int j = 0; j < 2; j++)
			{
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

	private void step()
	{
		long endStep = 0;
		long startStep = System.currentTimeMillis();
		long end = 0;
		long start = System.currentTimeMillis();

		logger.info("SPH clear buffer");
		runClearBuffers();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.CLEAR_BUFFERS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH clear buffer end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH hash particles");
		CLEvent hashParticles = runHashParticles();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.HASH_PARTICLES);
		}
		end = System.currentTimeMillis();
		logger.info("SPH hash particles end, took " + (end - start) + "ms");
		start = end;

		// host needs to wait as the next operation requires values from buffers
		hashParticles.waitFor();

		logger.info("SPH sort");
		runSort();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.SORT);
		}
		end = System.currentTimeMillis();
		logger.info("SPH sort end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH sort post pass");
		runSortPostPass();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.SORT_POST_PASS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH sort post pass end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH index");
		CLEvent runIndexx = runIndexx();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.INDEX);
		}
		end = System.currentTimeMillis();
		logger.info("SPH index end, took " + (end - start) + "ms");
		start = end;

		// host needs to wait as the next operation requires values from buffers
		runIndexx.waitFor();

		logger.info("SPH index post pass");
		runIndexPostPass();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.INDEX_POST_PASS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH index post pass end, took " + (end - start) + "ms");
		start = end;

		logger.info("SPH find neighbors");
		runFindNeighbors();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.FIND_NEIGHBORS);
		}
		end = System.currentTimeMillis();
		logger.info("SPH find neighbors end, took " + (end - start) + "ms");
		start = end;

		// PCISPH stuff starts here
		logger.info("PCI-SPH compute density");
		run_pcisph_computeDensity();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.COMPUTE_DENSITY);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH compute density end, took " + (end - start) + "ms");
		start = end;

		logger.info("PCI-SPH compute forces and init pressure");
		run_pcisph_computeForcesAndInitPressure();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.COMPUTE_FORCES_INIT_PRESSURE);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH compute forces and init pressure end, took " + (end - start) + "ms");
		start = end;

		// Do elastic stuff only if we have elastic particles
		if(_numOfElasticP > 0)
		{
			logger.info("PCI-SPH compute elastic forces");
			run_pcisph_computeElasticForces();
			if(_recordCheckPoints)
			{
				recordCheckpoints(KernelsEnum.COMPUTE_ELASTIC_FORCES);
			}
			end = System.currentTimeMillis();
			logger.info("PCI-SPH compute elastic forces end, took " + (end - start) + "ms");
			start = end;
		}

		logger.info("PCI-SPH predict/correct loop");
		// LOOP: 3 times or until "error" becomes less than 2%
		int iter = 0;
		int maxIterations = 3;
		do
		{
			run_pcisph_predictPositions();
			run_pcisph_predictDensity();
			run_pcisph_correctPressure();
			run_pcisph_computePressureForceAcceleration();

			iter++;
		}
		while((iter < maxIterations));
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.PREDICTIVE_LOOP);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH predict/correct loop end, took " + (end - start) + "ms");
		start = end;

		logger.info("PCI-SPH integrate");
		CLEvent event = run_pcisph_integrate();
		if(_recordCheckPoints)
		{
			recordCheckpoints(KernelsEnum.INTEGRATE);
		}
		end = System.currentTimeMillis();
		logger.info("PCI-SPH integrate end, took " + (end - start) + "ms");
		start = end;

		// wait for the end of the run_pcisph_integrate on device
		event.waitFor();

		logger.info("SPH finish queue");
		// TODO: figure out if we need to actually call this
		_queue.finish();
		end = System.currentTimeMillis();
		logger.info("SPH finish queue end, took " + (end - start) + "ms");
		start = end;

		endStep = System.currentTimeMillis();
		logger.info("SPH computation step done, took " + (endStep - startStep) + "ms");
	}

	public void finishQueue()
	{
		_queue.finish();
	}

	private int getParticleCountRoundedUp()
	{
		return (((_particleCount - 1) / 256) + 1) * 256;
	}

	@Override
	public void solve(IAspectConfiguration aspectConfiguration, AspectNode aspect) throws GeppettoExecutionException
	{
		// TODO: extend this to use time configuration to do multiple steps in
		// one go
		long time = System.currentTimeMillis();
		logger.info("SPH solver start");

		// TODO: setActivationSignal

		long end = 0;
		long start = System.currentTimeMillis();
		logger.info("SPH STEP START");
		step();
		updateStateTree(aspect);

		end = System.currentTimeMillis();
		logger.info("SPH STEP END, took " + (end - start) + "ms");

		logger.info("SPH solver end, took: " + (System.currentTimeMillis() - time) + "ms");
	}

	private void updateStateTree(AspectNode aspect)
	{
		AspectSubTreeNode simulationTree = (AspectSubTreeNode) aspect.getSubTree(AspectTreeType.SIMULATION_TREE);

		_positionPtr = _position.map(_queue, CLMem.MapFlags.Read);

		// ASSUMPTION: The solver will never create new states after the first
		// time step
		// we can call it principle of conservation of the states; if there is a
		// good
		// reason to revoke this assumption we need to add code that at every
		// cycle checks
		// if some new states exist to eventually add them to the stateTree

		updateSimulationTree(simulationTree);

		_position.unmap(_queue, _positionPtr);

		simulationTree.setModified(true);
		AspectNode aspectNode = (AspectNode) simulationTree.getParent();
		aspectNode.setModified(true);
		if(aspectNode.getParent() != null)
		{
			((EntityNode) aspectNode.getParentEntity()).updateParentEntitiesFlags(true);
		}

	}

	/**
	 * Updates nodes for simulation tree, if tree is empty call method to populate it.
	 * 
	 * @param simulationTree
	 */
	private void updateSimulationTree(AspectSubTreeNode simulationTree)
	{
		if(simulationTree.getChildren().isEmpty())
		{
			simulationTree.setId(AspectTreeType.SIMULATION_TREE.toString());
		}
		else
		{
			// watch tree not empty populate new values
			UpdateSPHSimulationTreeVisitor visitor = new UpdateSPHSimulationTreeVisitor(_positionPtr);
			simulationTree.apply(visitor);
		}
	}

	private boolean containsNode(ACompositeValue node, String name)
	{
		List<ANode> children = node.getChildren();

		boolean addNewNode = true;
		for(ANode child : children)
		{
			if(child.getId().equals(name))
			{
				addNewNode = false;
				return addNewNode;
			}
			if(child instanceof ACompositeValue)
			{
				if(((ACompositeValue) child).getChildren() != null)
				{
					addNewNode = containsNode((ACompositeValue) child, name);
				}
			}

		}

		return addNewNode;
	}

	private ACompositeValue getNode(ACompositeValue node, String name)
	{
		ACompositeValue newNode = null;

		List<ANode> children = node.getChildren();
		for(ANode child : children)
		{
			if(child.getId().equals(name))
			{
				newNode = (ACompositeValue) child;
				return newNode;
			}
			if(child instanceof ACompositeValue)
			{
				if(((ACompositeValue) child).getChildren() != null)
				{
					newNode = getNode((ACompositeValue) child, name);
				}
			}

		}

		return newNode;
	}

	@Override
	public void initialize(IModel model)
	{
		_model = (SPHModelX) model;

		setBuffersFromModel();
	}

	@Override
	public void dispose()
	{
		// close the context and "buonanotte al secchio" (good night to the
		// bucket)
		cleanContext();
	}

	private void setActivationSignal(float[] activation)
	{
		// put results back
		_activationSignalPtr = _activationSignal.map(_queue, CLMem.MapFlags.Write);
		_activationSignalPtr.setFloats(activation);
		_activationSignal.unmap(_queue, _activationSignalPtr);
	}

	private void recordCheckpoints(KernelsEnum kernelCheckpoint)
	{
		PCISPHCheckPoint check = new PCISPHCheckPoint();

		// read buffers into lists and populate checkpoint object
		check.acceleration = this.<Float> getBufferValues(_accelerationPtr, _acceleration, this._buffersSizeMap.get(BuffersEnum.ACCELERATION));
		check.gridCellIndex = this.<Integer> getBufferValues(_gridCellIndexPtr, _gridCellIndex, this._buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX));
		check.gridCellIndexFixedUp = this.<Integer> getBufferValues(_gridCellIndexFixedUpPtr, _gridCellIndexFixedUp, this._buffersSizeMap.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		check.neighborMap = this.<Float> getBufferValues(_neighborMapPtr, _neighborMap, this._buffersSizeMap.get(BuffersEnum.NEIGHBOR_MAP));
		check.particleIndex = this.<Integer> getBufferValues(_particleIndexPtr, _particleIndex, this._buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX));
		check.particleIndexBack = this.<Integer> getBufferValues(_particleIndexBackPtr, _particleIndexBack, this._buffersSizeMap.get(BuffersEnum.PARTICLE_INDEX_BACK));
		check.position = this.<Float> getBufferValues(_positionPtr, _position, this._buffersSizeMap.get(BuffersEnum.POSITION));
		check.pressure = this.<Float> getBufferValues(_pressurePtr, _pressure, this._buffersSizeMap.get(BuffersEnum.PRESSURE));
		check.rho = this.<Float> getBufferValues(_rhoPtr, _rho, this._buffersSizeMap.get(BuffersEnum.RHO));
		check.sortedPosition = this.<Float> getBufferValues(_sortedPositionPtr, _sortedPosition, this._buffersSizeMap.get(BuffersEnum.SORTED_POSITION));
		check.sortedVelocity = this.<Float> getBufferValues(_sortedVelocityPtr, _sortedVelocity, this._buffersSizeMap.get(BuffersEnum.SORTED_VELOCITY));
		check.velocity = this.<Float> getBufferValues(_velocityPtr, _velocity, this._buffersSizeMap.get(BuffersEnum.VELOCITY));
		if(_numOfElasticP > 0)
		{
			check.elasticConnections = this.<Float> getBufferValues(_elasticConnectionsDataPtr, _elasticConnectionsData, this._buffersSizeMap.get(BuffersEnum.ELASTIC_CONNECTIONS));
		}

		_checkpointsMap.put(kernelCheckpoint, check);
	}

	/*
	 * A method to retrieve buffer values into simple lists
	 */
	private <T> List<T> getBufferValues(Pointer<T> pointer, CLBuffer<T> buffer, int size)
	{
		List<T> list = new ArrayList<T>();

		pointer = buffer.map(_queue, CLMem.MapFlags.Read);

		for(int i = 0; i < size; i++)
		{
			list.add(pointer.get(i));
		}

		buffer.unmap(_queue, pointer);

		return list;
	}

	@Override
	public void updateVisualizationTree(AspectNode aspect)
	{
		AspectSubTreeNode visualTree = (AspectSubTreeNode) aspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);

		UpdateSPHVisualizationTreeVisitor updateSPHStateTreeVisitor = new UpdateSPHVisualizationTreeVisitor(_positionPtr);
		visualTree.apply(updateSPHStateTreeVisitor);

		visualTree.setModified(true);
	}

	@Override
	public void registerGeppettoService() throws Exception
	{

	}

};
