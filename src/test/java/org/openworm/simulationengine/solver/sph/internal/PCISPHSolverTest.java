/**
 * 
 */
package org.openworm.simulationengine.solver.sph.internal;

import static org.junit.Assert.fail;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.junit.Test;
import org.openworm.simulationengine.core.model.IModel;
import org.openworm.simulationengine.model.sph.services.SPHModelInterpreterService;
import org.openworm.simulationengine.model.sph.x.SPHModelX;
import org.openworm.simulationengine.solver.sph.SPHSolverService;

/**
 * @author matteo@openworm.org
 * @author giovanni@openworm.org
 */
public class PCISPHSolverTest
{
	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoCrash1()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/eshuozw196k3vci/sphModel_14.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 1; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoCrash2()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/8869zlz971ogyra/sphModel_small.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 18; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 * Same scene as testSolve14 but with 1 more particle
	 */
	@Test
	public void testSolve15_NoCrash()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/9kx2p8qspdgphd4/sphModel_15.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 254; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 * 296 boundary particles + 216 liquid particles (total of 512 particles) with random position and velocity
	 */
	@Test
	public void testSolve216_NoCrash()
	{
		try
		{
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			URL url = new URL("https://www.dropbox.com/s/lerz4rkx75nq0bk/sphModel_216.xml?dl=1");
			List<IModel> models = modelInterpreter.readModel(url);
			for (int cycles = 0; cycles < 1; cycles++)
			{
				models=solver.solve(models, null).get(0);
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 *  A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolvePureLiquidScene_ParticlesMoving()
	{
		try
		{
			URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");
			
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			
			List<IModel> initial_models = modelInterpreter.readModel(url);
			List<IModel> models =  new ArrayList<IModel>(initial_models);
			
			for (int cycles = 0; cycles < 1; cycles++)
			{
				models=solver.solve(models, null).get(0);
				int pd=((SPHModelX)initial_models.get(0)).compareTo((SPHModelX)(models.get(0)));
				System.out.println("Particles different at cycle "+cycles+": "+pd);
				Assert.assertFalse(pd==0);
				initial_models=models;
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 *  A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolveElastic_NoCrash()
	{
		try
		{
			URL url = this.getClass().getResource("/sphModel_Elastic.xml");
			
			SPHSolverService solver = new SPHSolverService();
			SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
			
			List<IModel> initial_models = modelInterpreter.readModel(url);
			List<IModel> models =  new ArrayList<IModel>(initial_models);
			
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
}
