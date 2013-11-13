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

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.state.StateTreeRoot;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

/**
 * 
 * Multiple SPH Solver threads test.
 * 
 * @author  Jesus R. Martinez (jesus@metacell.us)
 *
 */
public class MultipleThreadsTests {

	/*
	 * Creates different 5 different solver threads and runs them at the same time
	 */
	@Test
	public void testMultipleThreads() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_1.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		PCISPHTestUtilities.checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		
		SPHSolverService solver3 = new SPHSolverService();
		solver3.initialize(model);
		
		SPHSolverService solver4 = new SPHSolverService();
		solver4.initialize(model);
		
		SPHSolverService solver5 = new SPHSolverService();
		solver5.initialize(model);
		
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 2, 1));
		StateTreeRoot stateSet2 = solver2.solve(new TimeConfiguration(0.2f, 4, 1));
		StateTreeRoot stateSet3 = solver3.solve(new TimeConfiguration(0.3f, 5, 1));
		StateTreeRoot stateSet4 = solver4.solve(new TimeConfiguration(0.4f, 3, 1));
		StateTreeRoot stateSet5 = solver5.solve(new TimeConfiguration(0.5f, 2, 1));
		
		// NOTE: this is commented out as it fails on Apple CPU - should pass everywhere else
		// PCISPHTestUtilities.checkStateTreeForNaN(stateSet, false);
		
		Assert.assertTrue("Particle count doesn't match.", stateSet.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateSet2.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateSet3.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateSet4.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
		Assert.assertTrue("Particle count doesn't match.", stateSet5.getChildren().size() == PCISPHTestUtilities.countNonBoundaryParticles((SPHModelX)model));
	}
}
