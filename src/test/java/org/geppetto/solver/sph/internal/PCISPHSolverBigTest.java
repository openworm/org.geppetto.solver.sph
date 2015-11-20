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

import org.geppetto.core.data.model.IAspectConfiguration;
import org.geppetto.core.data.model.local.LocalAspectConfiguration;
import org.geppetto.core.data.model.local.LocalInstancePath;
import org.geppetto.core.data.model.local.LocalSimulatorConfiguration;
import org.geppetto.core.model.IModel;
import org.geppetto.core.model.runtime.AspectSubTreeNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode.AspectTreeType;
import org.geppetto.core.model.typesystem.AspectNode;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class PCISPHSolverBigTest {
	
	LocalSimulatorConfiguration sc=new LocalSimulatorConfiguration(1, "", "", 0.0002f, 0.2f, null);
	LocalInstancePath ip=new LocalInstancePath(1, "", "", "");
	IAspectConfiguration aspectConfiguration=new LocalAspectConfiguration(1, ip, null, null, sc);
	
	/*
	 * A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolvePureLiquidScene_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// NOTE: commenting out as it takes a while for big scenes
		//checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		AspectNode stateSet = new AspectNode("StateSet");
		solver.solve(aspectConfiguration,stateSet);
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateSet.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateSet.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
	
	/*
	 * Pure liquid scene with "many" liquid particles
	 * NOTE: compares results from running cycles in one go vs step by step
	 */
	@Test
	public void testSolvePureLiquid_StepByStep_VS_OneGo() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		//checkModelForOverlappingParticles((SPHModelX)model, false);
		
		int cycles = 2;
		
		// run cycles one by one
		SPHSolverService solver1 = new SPHSolverService();
		solver1.initialize(model);
		AspectNode stateTree1 = new AspectNode("StateSet1");
		for(int i = 0; i < cycles; i++){
			solver1.solve(aspectConfiguration,stateTree1);
		}
		
		// run cycles at once
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		
		AspectNode stateTree2 = new AspectNode("StateSet2");
		solver2.solve(aspectConfiguration,stateTree2);
		
		//checks the trees are equivalent
		Assert.assertEquals(stateTree1.toString(),stateTree2.toString());
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateTree1.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateTree2.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateTree1.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateTree2.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}

	/*
	 * Scene with a block of elastic matter plus liquid
	 */
	@Test
	public void testSolveElastic_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		//checkModelForOverlappingParticles((SPHModelX)model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		AspectNode stateSet = new AspectNode("StateSet");
		solver.solve(aspectConfiguration,stateSet);
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateSet.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateSet.getSubTree(AspectSubTreeNode.AspectTreeType.MODEL_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
	
	/*
	 * Scene with a block of elastic matter plus liquid
	 * NOTE: solves step by step - this frees memory so we can run many steps without running out of heap space
	 */
	@Test
	public void testSolveElastic_StepByStep_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		//checkModelForOverlappingParticles((SPHModelX)model, false);
		
		int cycles = 2;
		
		// run cycles one by one
		SPHSolverService solver1 = new SPHSolverService();
		solver1.initialize(model);		
		AspectNode stateTree = new AspectNode("State Tree");
		for(int i = 0; i < cycles; i++){
			solver1.solve(aspectConfiguration,stateTree);
			if (i > 0) {
				PCISPHTestUtilities.removeFirstStateFromTree((AspectSubTreeNode) stateTree.getSubTree(AspectTreeType.VISUALIZATION_TREE));
			}
		}
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateTree.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateTree.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
	
	/*
	 * Scene with a block of elastic matter plus liquid
	 * NOTE: compares results from running cycles in one go vs step by step
	 */
	@Test
	public void testSolveElastic_StepByStep_VS_OneGo() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");
		
		// check that we don't have particles with overlapping positions
		//checkModelForOverlappingParticles((SPHModelX)model, false);
		
		int cycles = 2;
		
		// run cycles one by one
		SPHSolverService solver1 = new SPHSolverService();
		solver1.initialize(model);
		AspectNode stateTree1 = new AspectNode("Aspect");
		for(int i = 0; i < cycles; i++){
			solver1.solve(aspectConfiguration,stateTree1);
		}
		
		// run cycles at once
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		
		AspectNode stateTree2 = new AspectNode("Aspect");
		solver2.solve(aspectConfiguration,stateTree2);
		
		//checks the trees are equivalent
		Assert.assertEquals(stateTree1.toString(),stateTree2.toString());
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateTree1.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateTree2.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateTree1.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateTree2.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
	
	/*
	 * Scene with a block of elastic matter plus liquid
	 */
	@Test
	public void testSolveElastic_Contractible_7220_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_elastic_contractible_7220.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url,null,"");

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		AspectNode stateSet = new AspectNode("Aspect");
		solver.solve(aspectConfiguration,stateSet);
		
		PCISPHTestUtilities.checkStateTreeForNaN((AspectSubTreeNode) stateSet.getSubTree(AspectTreeType.VISUALIZATION_TREE), false);
		
		Assert.assertTrue("Particle count doesn't match.", stateSet.getSubTree(AspectSubTreeNode.AspectTreeType.VISUALIZATION_TREE).getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
}
