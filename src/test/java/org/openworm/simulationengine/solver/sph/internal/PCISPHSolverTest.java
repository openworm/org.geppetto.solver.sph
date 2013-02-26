/**
 * 
 */
package org.openworm.simulationengine.solver.sph.internal;

import static org.junit.Assert.*;

import java.net.URL;
import java.util.List;

import org.junit.Test;
import org.openworm.simulationengine.core.model.IModel;
import org.openworm.simulationengine.core.simulation.ITimeConfiguration;
import org.openworm.simulationengine.core.simulation.TimeConfiguration;
import org.openworm.simulationengine.model.sph.services.SPHModelInterpreterService;
import org.openworm.simulationengine.solver.sph.SPHSolverService;

/**
 * @author matteocantarelli
 * 
 */
public class PCISPHSolverTest
{

	/**
	 * Test method for {@link org.openworm.simulationengine.solver.sph.SPHSolverService#solve(java.util.List, org.openworm.simulationengine.core.simulation.ITimeConfiguration)}.
	 */
	@Test
	public void testSolve14()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/8869zlz971ogyra/sphModel_small.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 10; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/**
	 * Test method for {@link org.openworm.simulationengine.solver.sph.SPHSolverService#solve(java.util.List, org.openworm.simulationengine.core.simulation.ITimeConfiguration)}.
	 */
	@Test
	public void testSolve1024()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/yvhca5hwx9uyw5o/sphModel_1024.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 23; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/**
	 * Test method for {@link org.openworm.simulationengine.solver.sph.SPHSolverService#solve(java.util.List, org.openworm.simulationengine.core.simulation.ITimeConfiguration)}.
	 */
	@Test
	public void testSolve2048()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/deb035iisfcdabp/sphModel_2048.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 23; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}

}
