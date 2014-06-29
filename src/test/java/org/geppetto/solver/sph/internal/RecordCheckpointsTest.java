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

package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.Map;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.state.AspectTreeNode;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.BuffersEnum;
import org.geppetto.solver.sph.KernelsEnum;
import org.geppetto.solver.sph.PCISPHCheckPoint;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class RecordCheckpointsTest {

	@Test
	public void testRecordCheckpointsIsEmpty() throws Exception {
		URL url = this.getClass().getResource("/sphModel_15.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX)model, false);
		
		// use default constructor - doesn't record checkpoints
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		AspectTreeNode stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
		PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
		Assert.assertTrue("Particle count doesn't match.", stateSet.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		
		Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		Assert.assertTrue("The map of checkpoints should be empty but is not", checkpoints.isEmpty());
	}
	
	@Test
	public void testRecordCheckpoints() throws Exception {
		URL url = this.getClass().getResource("/sphModel_15.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX)model, false);
		
		// tell the solver to record checkpoints
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		AspectTreeNode stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
		PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
		Assert.assertTrue("Particle count doesn't match.", stateSet.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		
		Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		
		Assert.assertFalse("The map of checkpoints should not be empty but it is", checkpoints.isEmpty());
		
		Map<BuffersEnum, Integer> bufferSize = solver.getBuffersSizeMap();
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint firstCheckPoint = checkpoints.get(KernelsEnum.CLEAR_BUFFERS);
		Assert.assertTrue("1st checkpoint: size of acceleration buffer doesn't match expected value", firstCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("1st checkpoint: size of gridCellIndex buffer doesn't match expected value", firstCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("1st checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", firstCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("1st checkpoint: size of neighborMap buffer doesn't match expected value", firstCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("1st checkpoint: size of particleIndex buffer doesn't match expected value", firstCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("1st checkpoint: size of particelIndexBack buffer doesn't match expected value", firstCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("1st checkpoint: size of position buffer doesn't match expected value", firstCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("1st checkpoint: size of pressure buffer doesn't match expected value", firstCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("1st checkpoint: size of rho buffer doesn't match expected value", firstCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("1st checkpoint: size of sortedPosition buffer doesn't match expected value", firstCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("1st checkpoint: size of sortedVelocity buffer doesn't match expected value", firstCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("1st checkpoint: size of velocity buffer doesn't match expected value", firstCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		// no elastic connections in the model loaded
		Assert.assertTrue("1st checkpoint: elastic connections buffer recordings should be empty but they're not", firstCheckPoint.elasticConnections.isEmpty());
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint lastCheckPoint = checkpoints.get(KernelsEnum.INTEGRATE);
		Assert.assertTrue("Last checkpoint: size of acceleration buffer doesn't match expected value", lastCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("Last checkpoint: size of gridCellIndex buffer doesn't match expected value", lastCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("Last checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", lastCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("Last checkpoint: size of neighborMap buffer doesn't match expected value", lastCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("Last checkpoint: size of particleIndex buffer doesn't match expected value", lastCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("Last checkpoint: size of particelIndexBack buffer doesn't match expected value", lastCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("Last checkpoint: size of position buffer doesn't match expected value", lastCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("Last checkpoint: size of pressure buffer doesn't match expected value", lastCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("Last checkpoint: size of rho buffer doesn't match expected value", lastCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("Last checkpoint: size of sortedPosition buffer doesn't match expected value", lastCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("Last checkpoint: size of sortedVelocity buffer doesn't match expected value", lastCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("Last checkpoint: size of velocity buffer doesn't match expected value", lastCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		// no elastic connections in the model loaded
		Assert.assertTrue("Last checkpoint: elastic connections buffer recordings should be empty but they're not", lastCheckPoint.elasticConnections.isEmpty());
	}
	
	@Test
	public void testRecordCheckpointsMultiStep() throws Exception {
		URL url = this.getClass().getResource("/sphModel_15.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX)model, false);
		
		// tell the solver to record checkpoints
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		// run many steps - test should give same results
		AspectTreeNode stateSet = solver.solve(new TimeConfiguration(0.1f, 10, 1));
		
		PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
		Assert.assertTrue("Particle count doesn't match.", stateSet.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		
		Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		
		Assert.assertFalse("The map of checkpoints should not be empty but it is", checkpoints.isEmpty());
		
		Map<BuffersEnum, Integer> bufferSize = solver.getBuffersSizeMap();
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint firstCheckPoint = checkpoints.get(KernelsEnum.CLEAR_BUFFERS);
		Assert.assertTrue("1st checkpoint: size of acceleration buffer doesn't match expected value", firstCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("1st checkpoint: size of gridCellIndex buffer doesn't match expected value", firstCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("1st checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", firstCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("1st checkpoint: size of neighborMap buffer doesn't match expected value", firstCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("1st checkpoint: size of particleIndex buffer doesn't match expected value", firstCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("1st checkpoint: size of particelIndexBack buffer doesn't match expected value", firstCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("1st checkpoint: size of position buffer doesn't match expected value", firstCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("1st checkpoint: size of pressure buffer doesn't match expected value", firstCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("1st checkpoint: size of rho buffer doesn't match expected value", firstCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("1st checkpoint: size of sortedPosition buffer doesn't match expected value", firstCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("1st checkpoint: size of sortedVelocity buffer doesn't match expected value", firstCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("1st checkpoint: size of velocity buffer doesn't match expected value", firstCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		// no elastic connections in the model loaded
		Assert.assertTrue("1st checkpoint: elastic connections buffer recordings should be empty but they're not", firstCheckPoint.elasticConnections.isEmpty());
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint lastCheckPoint = checkpoints.get(KernelsEnum.INTEGRATE);
		Assert.assertTrue("Last checkpoint: size of acceleration buffer doesn't match expected value", lastCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("Last checkpoint: size of gridCellIndex buffer doesn't match expected value", lastCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("Last checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", lastCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("Last checkpoint: size of neighborMap buffer doesn't match expected value", lastCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("Last checkpoint: size of particleIndex buffer doesn't match expected value", lastCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("Last checkpoint: size of particelIndexBack buffer doesn't match expected value", lastCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("Last checkpoint: size of position buffer doesn't match expected value", lastCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("Last checkpoint: size of pressure buffer doesn't match expected value", lastCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("Last checkpoint: size of rho buffer doesn't match expected value", lastCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("Last checkpoint: size of sortedPosition buffer doesn't match expected value", lastCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("Last checkpoint: size of sortedVelocity buffer doesn't match expected value", lastCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("Last checkpoint: size of velocity buffer doesn't match expected value", lastCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		// no elastic connections in the model loaded
		Assert.assertTrue("Last checkpoint: elastic connections buffer recordings should be empty but they're not", lastCheckPoint.elasticConnections.isEmpty());
	}
	
	@Test
	public void testRecordCheckpointsWithElasticConnections() throws Exception {
		URL url = this.getClass().getResource("/sphModel_elastic_contractible_7220.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
	
		// tell the solver to record checkpoints
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		AspectTreeNode stateSet = solver.solve(new TimeConfiguration(null, 1, null));
		
		PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
		
		Assert.assertTrue("Particle count doesn't match.", stateSet.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		
		Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		Assert.assertFalse("The map of checkpoints should not be empty but it is", checkpoints.isEmpty());
		
		Map<BuffersEnum, Integer> bufferSize = solver.getBuffersSizeMap();
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint firstCheckPoint = checkpoints.get(KernelsEnum.CLEAR_BUFFERS);
		Assert.assertTrue("1st checkpoint: size of acceleration buffer doesn't match expected value", firstCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("1st checkpoint: size of gridCellIndex buffer doesn't match expected value", firstCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("1st checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", firstCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("1st checkpoint: size of neighborMap buffer doesn't match expected value", firstCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("1st checkpoint: size of particleIndex buffer doesn't match expected value", firstCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("1st checkpoint: size of particelIndexBack buffer doesn't match expected value", firstCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("1st checkpoint: size of position buffer doesn't match expected value", firstCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("1st checkpoint: size of pressure buffer doesn't match expected value", firstCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("1st checkpoint: size of rho buffer doesn't match expected value", firstCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("1st checkpoint: size of sortedPosition buffer doesn't match expected value", firstCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("1st checkpoint: size of sortedVelocity buffer doesn't match expected value", firstCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("1st checkpoint: size of velocity buffer doesn't match expected value", firstCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		Assert.assertTrue("1st checkpoint: size of elasticConnections buffer doesn't match expected value", firstCheckPoint.elasticConnections.size() == bufferSize.get(BuffersEnum.ELASTIC_CONNECTIONS));
		
		// check buffer sizes at first and last checkpoints
		PCISPHCheckPoint lastCheckPoint = checkpoints.get(KernelsEnum.INTEGRATE);
		Assert.assertTrue("Last checkpoint: size of acceleration buffer doesn't match expected value", lastCheckPoint.acceleration.size() == bufferSize.get(BuffersEnum.ACCELERATION));
		Assert.assertTrue("Last checkpoint: size of gridCellIndex buffer doesn't match expected value", lastCheckPoint.gridCellIndex.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX));
		Assert.assertTrue("Last checkpoint: size of gridCellIndexFixedUp buffer doesn't match expected value", lastCheckPoint.gridCellIndexFixedUp.size() == bufferSize.get(BuffersEnum.GRID_CELL_INDEX_FIXED));
		Assert.assertTrue("Last checkpoint: size of neighborMap buffer doesn't match expected value", lastCheckPoint.neighborMap.size() == bufferSize.get(BuffersEnum.NEIGHBOR_MAP));
		Assert.assertTrue("Last checkpoint: size of particleIndex buffer doesn't match expected value", lastCheckPoint.particleIndex.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX));
		Assert.assertTrue("Last checkpoint: size of particelIndexBack buffer doesn't match expected value", lastCheckPoint.particleIndexBack.size() == bufferSize.get(BuffersEnum.PARTICLE_INDEX_BACK));
		Assert.assertTrue("Last checkpoint: size of position buffer doesn't match expected value", lastCheckPoint.position.size() == bufferSize.get(BuffersEnum.POSITION));
		Assert.assertTrue("Last checkpoint: size of pressure buffer doesn't match expected value", lastCheckPoint.pressure.size() == bufferSize.get(BuffersEnum.PRESSURE));
		Assert.assertTrue("Last checkpoint: size of rho buffer doesn't match expected value", lastCheckPoint.rho.size() == bufferSize.get(BuffersEnum.RHO));
		Assert.assertTrue("Last checkpoint: size of sortedPosition buffer doesn't match expected value", lastCheckPoint.sortedPosition.size() == bufferSize.get(BuffersEnum.SORTED_POSITION));
		Assert.assertTrue("Last checkpoint: size of sortedVelocity buffer doesn't match expected value", lastCheckPoint.sortedVelocity.size() == bufferSize.get(BuffersEnum.SORTED_VELOCITY));
		Assert.assertTrue("Last checkpoint: size of velocity buffer doesn't match expected value", lastCheckPoint.velocity.size() == bufferSize.get(BuffersEnum.VELOCITY));
		Assert.assertTrue("Last checkpoint: size of elasticConnections buffer doesn't match expected value", lastCheckPoint.elasticConnections.size() == bufferSize.get(BuffersEnum.ELASTIC_CONNECTIONS));
	}
}
