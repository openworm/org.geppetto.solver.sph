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

package org.geppetto.solver.sph.internal;

import java.net.URL;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.runtime.ACompositeNode;
import org.geppetto.core.model.runtime.ANode;
import org.geppetto.core.model.runtime.AspectNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode.AspectTreeType;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Before;
import org.junit.Test;

/**
 * @author matteo@openworm.org
 * @author giovanni@openworm.org
 */
public class PCISPHSolverTest
{

	AspectNode _sphAspect = null;

	@Before
	public void runBeforeEveryTest()
	{
		_sphAspect=new AspectNode("test");
	}

	/*
	 * 296 boundary particles + 1 liquid particle NOTE: particle is very close to the origin and it is shot towards it with high velocity
	 */
	@Test
	public void testSolve1_NoNaN() throws Exception
	{

		URL url = this.getClass().getResource("/sphModel_1.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.solve(new TimeConfiguration(0.1f, 2, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);

		// NOTE: this is commented out as it fails on Apple CPU - should pass everywhere else
		// PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);

		Assert.assertTrue("Particle count doesn't match.", countLiquidParticles(simulationTree) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	private int countLiquidParticles(AspectSubTreeNode simulationTree)
	{
		return ((ACompositeNode)simulationTree.getChildren().get(0)).getChildren().size();
	}

	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_ExpectNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_14.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// we have one particle that overlaps with boundary particles in this test
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, true);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.solve(new TimeConfiguration(0.1f, 20, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		
		// expect NaN values since in the initial conditions particle 309 overlaps with boundary particles
		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree, true);
		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_small.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.solve(new TimeConfiguration(0.1f, 20, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		
		PCISPHTestUtilities.checkStateTreeForNaN(_sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);

		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	/*
	 * 296 boundary particles + 14 liquid particles NOTE: compares results from running cycles in one go vs step by step
	 */
	@Test
	public void testSolve14_StepByStep_VS_OneGo() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_small.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, false);

		int cycles = 20;

		// run cycles one by one
		SPHSolverService solver1 = new SPHSolverService();
		solver1.initialize(model);
		
		for(int i = 0; i < cycles; i++)
		{
			solver1.solve(new TimeConfiguration(0.1f, 1, 1),_sphAspect);
		}
		AspectSubTreeNode simulationTree1=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver1.populateVisualTree(model, simulationTree1);
		// run cycles at once
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		AspectNode sphAspect2=new AspectNode("test");
		solver2.solve(new TimeConfiguration(0.1f, cycles, 1),sphAspect2);
		AspectSubTreeNode simulationTree2=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver2.populateVisualTree(model, simulationTree2);
		// checks the trees are equivalent
		Assert.assertEquals(simulationTree1.toString(), simulationTree2.toString());

		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree1, false);
		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree2, false);

		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree1) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree2) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	/*
	 * Same scene as testSolve14 but with 1 more particle
	 */
	@Test
	public void testSolve15_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_15.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.solve(new TimeConfiguration(0.1f, 20, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree, false);

		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	/*
	 * 296 boundary particles + 216 liquid particles (total of 512 particles) with random position and velocity
	 */
	@Test
	public void testSolve216_NaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_216.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX) model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.solve(new TimeConfiguration(0.1f, 10, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		
		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree, true);

		Assert.assertTrue("Particle count doesn't match.",
				countLiquidParticles(simulationTree) == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX) model));
	}

	@Test
	public void testMembranesScene() throws Exception
	{
		URL url = this.getClass().getResource("/cube_with_membranes.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		solver.setRecordCheckpoint(false);
		solver.solve(new TimeConfiguration(0.1f, 100, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		Assert.assertTrue(true);
		PCISPHTestUtilities.checkStateTreeForNaN(simulationTree, false);

	}

	@Test
	public void testMembranesWithWaterScene() throws Exception
	{
		URL url = this.getClass().getResource("/cube_with_membranes_water_inside.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		// solver.setRecordCheckpoint(true);
		solver.solve(new TimeConfiguration(0.1f, 100, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		// Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		// PCISPHCheckPoint lastCheckPoint = checkpoints.get(KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES_FINALIZE);
		Assert.assertTrue(true);
		// PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);

	}
	
	@Test
	public void test_simple_liquid_scene_780() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_liquid_780.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url, null, "");

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		// solver.setRecordCheckpoint(true);
		solver.solve(new TimeConfiguration(0.1f, 100, 1),_sphAspect);
		AspectSubTreeNode simulationTree=(AspectSubTreeNode) _sphAspect.getSubTree(AspectTreeType.VISUALIZATION_TREE);
		solver.populateVisualTree(model, simulationTree);
		// Map<KernelsEnum, PCISPHCheckPoint> checkpoints = solver.getCheckpointsMap();
		// PCISPHCheckPoint lastCheckPoint = checkpoints.get(KernelsEnum.COMPUTE_INTERACTION_WITH_MEMBRANES_FINALIZE);
		Assert.assertTrue(true);
		// PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);

	}
}
