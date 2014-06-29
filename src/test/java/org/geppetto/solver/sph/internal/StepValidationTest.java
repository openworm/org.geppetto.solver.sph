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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import junit.framework.Assert;

import org.geppetto.core.model.state.AspectTreeNode;
import org.geppetto.core.simulation.TimeConfiguration;

import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class StepValidationTest {

	@Test
	public void testValidateLiquidScene16974() throws Exception {
		// load reference values at various steps from C++ version
		String position0 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_0.txt").getPath());
		String position1 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_1.txt").getPath());
		String position2 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_2.txt").getPath());
		String position3 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_3.txt").getPath());
		String position4 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_4.txt").getPath());
		String position5 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_5.txt").getPath());
		String position10 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_10.txt").getPath());
		String position20 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_20.txt").getPath());
		String position30 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_30.txt").getPath());
		String position40 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_40.txt").getPath());
		String position50 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_50.txt").getPath());
		
		Map<Integer, String[]> referenceValuesMap = new HashMap<Integer, String[]>();
		referenceValuesMap.put(0, position0.split(System.getProperty("line.separator")));
		referenceValuesMap.put(1, position1.split(System.getProperty("line.separator")));
		referenceValuesMap.put(2, position2.split(System.getProperty("line.separator")));
		referenceValuesMap.put(3, position3.split(System.getProperty("line.separator")));
		referenceValuesMap.put(4, position4.split(System.getProperty("line.separator")));
		referenceValuesMap.put(5, position5.split(System.getProperty("line.separator")));
		referenceValuesMap.put(10, position10.split(System.getProperty("line.separator")));
		referenceValuesMap.put(20, position20.split(System.getProperty("line.separator")));
		referenceValuesMap.put(30, position30.split(System.getProperty("line.separator")));
		referenceValuesMap.put(40, position40.split(System.getProperty("line.separator")));
		referenceValuesMap.put(50, position50.split(System.getProperty("line.separator")));
		
		// load Java generated scene
		URL url = this.getClass().getResource("/sphModel_liquid_16974.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		// 3. assert number of particles is fine
		for (Map.Entry<Integer, String[]> entry : referenceValuesMap.entrySet())
		{
		    Assert.assertTrue("number of lines on positions and number of particles on sphModel do not match", model.getParticles().size() == entry.getValue().length);
		}
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		Map<Integer,Set<Integer>> mismatchingSetsMap = new LinkedHashMap<Integer, Set<Integer>>();
		
		int step = 0;
		if( referenceValuesMap.containsKey(step) )
		{
			// get reference values
			String[] referenceValues = referenceValuesMap.get(step);

			AspectTreeNode stateSet = solver.getStateTree();

			// compare state tree with logged values for each recorded step
			CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
			stateSet.apply(compareVisitor);
			mismatchingSetsMap.put(step, compareVisitor.getMismatches());
		}
		
		for(int i = 0; i < 50; i++)
		{
			// calculate step
			AspectTreeNode stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
			
			// only keep latest step results
			if (i > 0) {
				PCISPHTestUtilities.removeFirstStateFromTree(stateSet);
			}
			
			// get reference values
			step = i + 1;
			String[] referenceValues = referenceValuesMap.get(step);
			
			PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
			
			if( referenceValuesMap.containsKey(step) )
			{
				// compare state tree with logged values for each recorded step
				CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
				stateSet.apply(compareVisitor);
				mismatchingSetsMap.put(step, compareVisitor.getMismatches());
			}
		}
		
		StringBuilder sb = new StringBuilder();
		
		// check mismatching stuff
		for (Map.Entry<Integer,Set<Integer>> entry : mismatchingSetsMap.entrySet())
		{
			Integer _step = entry.getKey();
			// assert there are no mismatching values
			if (entry.getValue().size() != 0)
				sb
					.append(entry.getValue().size()).append(" of ").append(referenceValuesMap.get(_step).length)
					.append(" mismatching values found on step ").append(_step).append("\n")
					;
		}
		
		String msg = sb.toString();
		if (!msg.isEmpty())
			Assert.fail(msg);
	}
	
	@Test
	public void testValidateElasticSceneSingleBundle() throws Exception {
		// load reference values at various steps from C++ version
		String position0 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/elastic/position_log_0.txt").getPath());
		String position1 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/elastic/position_log_1.txt").getPath());
		String position5 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/elastic/position_log_5.txt").getPath());
		String position10 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/elastic/position_log_10.txt").getPath());
		
		Map<Integer, String[]> referenceValuesMap = new HashMap<Integer, String[]>();
		referenceValuesMap.put(0, position0.split(System.getProperty("line.separator")));
		referenceValuesMap.put(1, position1.split(System.getProperty("line.separator")));
		referenceValuesMap.put(5, position5.split(System.getProperty("line.separator")));
		referenceValuesMap.put(10, position10.split(System.getProperty("line.separator")));
		
		// load Java generated scene
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		// 3. assert number of particles is fine
		for (Map.Entry<Integer, String[]> entry : referenceValuesMap.entrySet())
		{
		    Assert.assertTrue("number of lines on positions and number of particles on sphModel do not match", model.getParticles().size() == entry.getValue().length);
		}
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		Map<Integer,Set<Integer>> mismatchingSetsMap = new LinkedHashMap<Integer, Set<Integer>>();

		int step = 0;
		if( referenceValuesMap.containsKey(step) )
		{
			// get reference values
			String[] referenceValues = referenceValuesMap.get(step);
			
			AspectTreeNode stateSet = solver.getStateTree();

			// compare state tree with logged values for each recorded step
			CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
			stateSet.apply(compareVisitor);
			mismatchingSetsMap.put(step, compareVisitor.getMismatches());
		}

		for(int i = 0; i < 10; i++)
		{
			// calculate step
			AspectTreeNode stateSet = solver.solve(new TimeConfiguration(null, 1, null));
			
			// only keep latest step results
			if (i > 0) {
				PCISPHTestUtilities.removeFirstStateFromTree(stateSet);
			}
			
			// get reference values
			step = i + 1;
			String[] referenceValues = referenceValuesMap.get(step);
			
			PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
			
			if( referenceValuesMap.containsKey(step) )
			{
				// compare state tree with logged values for each recorded step
				CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
				stateSet.apply(compareVisitor);
				mismatchingSetsMap.put(step, compareVisitor.getMismatches());
			}
		}
		
		StringBuilder sb = new StringBuilder();
		
		// check mismatching stuff
		for (Map.Entry<Integer,Set<Integer>> entry : mismatchingSetsMap.entrySet())
		{
			Integer _step = entry.getKey();
			// assert there are no mismatching values
			if (entry.getValue().size() != 0)
				sb.append(entry.getValue().size()).append(" of ").append(referenceValuesMap.get(_step).length)
				  .append(" mismatching values found on step ").append(_step).append("\n");
		}
		
		String msg = sb.toString();
		if (!msg.isEmpty())
			Assert.fail(msg);
	}
	
	@Test
	public void testValidateLiquidScene780() throws Exception {
		// load reference values at various steps from C++ version
		String position0 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_0.txt").getPath());
		String position1 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_1.txt").getPath());
		String position2 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_2.txt").getPath());
		String position3 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_3.txt").getPath());
		String position4 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_4.txt").getPath());
		String position5 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_5.txt").getPath());
		String position10 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_10.txt").getPath());
		String position20 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_20.txt").getPath());
		String position30 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_30.txt").getPath());
		String position40 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_40.txt").getPath());
		String position50 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_50.txt").getPath());
		String position100 = PCISPHTestUtilities.readFile(StepValidationTest.class.getResource("/results/liquid_780/position_log_100.txt").getPath());
		
		Map<Integer, String[]> referenceValuesMap = new HashMap<Integer, String[]>();
		referenceValuesMap.put(0, position0.split(System.getProperty("line.separator")));
		referenceValuesMap.put(1, position1.split(System.getProperty("line.separator")));
		referenceValuesMap.put(2, position2.split(System.getProperty("line.separator")));
		referenceValuesMap.put(3, position3.split(System.getProperty("line.separator")));
		referenceValuesMap.put(4, position4.split(System.getProperty("line.separator")));
		referenceValuesMap.put(5, position5.split(System.getProperty("line.separator")));
		referenceValuesMap.put(10, position10.split(System.getProperty("line.separator")));
		referenceValuesMap.put(20, position20.split(System.getProperty("line.separator")));
		referenceValuesMap.put(30, position30.split(System.getProperty("line.separator")));
		referenceValuesMap.put(40, position40.split(System.getProperty("line.separator")));
		referenceValuesMap.put(50, position50.split(System.getProperty("line.separator")));
		referenceValuesMap.put(100, position100.split(System.getProperty("line.separator")));
		
		// load Java generated scene
		URL url = this.getClass().getResource("/sphModel_liquid_780.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		// 3. assert number of particles is fine
		for (Map.Entry<Integer, String[]> entry : referenceValuesMap.entrySet())
		{
		    Assert.assertTrue("number of lines on positions and number of particles on sphModel do not match", model.getParticles().size() == entry.getValue().length);
		}
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		Map<Integer,Set<Integer>> mismatchingSetsMap = new LinkedHashMap<Integer, Set<Integer>>();
		
		int step = 0;
		if( referenceValuesMap.containsKey(step) )
		{
			// get reference values
			String[] referenceValues = referenceValuesMap.get(step);

			AspectTreeNode stateSet = solver.getStateTree();

			// compare state tree with logged values for each recorded step
			CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
			stateSet.apply(compareVisitor);
			mismatchingSetsMap.put(step, compareVisitor.getMismatches());
		}
		
		for(int i = 0; i < 100; i++)
		{
			// calculate step
			AspectTreeNode stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
			
			// only keep latest step results
			if (i > 0) {
				PCISPHTestUtilities.removeFirstStateFromTree(stateSet);
			}
			
			// get reference values
			step = i + 1;
			String[] referenceValues = referenceValuesMap.get(step);
			
			PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
			
			if( referenceValuesMap.containsKey(step) )
			{
				// compare state tree with logged values for each recorded step
				CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
				stateSet.apply(compareVisitor);
				mismatchingSetsMap.put(step, compareVisitor.getMismatches());
			}
		}
		
		StringBuilder sb = new StringBuilder();
		
		// check mismatching stuff
		for (Map.Entry<Integer,Set<Integer>> entry : mismatchingSetsMap.entrySet())
		{
			Integer _step = entry.getKey();
			// assert there are no mismatching values
			if (entry.getValue().size() != 0)
				sb
					.append(entry.getValue().size()).append(" of ").append(referenceValuesMap.get(_step).length)
					.append(" mismatching values found on step ").append(_step).append("\n")
					;
		}
		
		String msg = sb.toString();
		if (!msg.isEmpty())
			Assert.fail(msg);
	}
}
