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
 * @author matteocantarelli
 * 
 */
public class PCISPHSolverTest
{
	/*
	 * 296 boundary particles + 14 liquid particles - Use the configuration from step 17 of the testSolve14 scene, crashes immediately
	 * NOTE: commented out not to break "maven install" build
	 */
	@Test
	public void testSolve14_ImmediateCrash()
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
	 * 296 boundary particles + 14 liquid particles - this runs fine for 18 steps then crashes
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
			for (int cycles = 0; cycles < 18; cycles++)
			{
				models=solver.solve(models, null).get(0);
				
				/*if(cycles == 17){
					
					for(SPHParticle p : ((SPHModel) models.get(0)).getParticles()){
						if(p.getPositionVector().getP() == SPHConstants.LIQUID_TYPE)
						System.out.println(p.getPositionVector().getX() + " / " + p.getPositionVector().getY() + " / " + p.getPositionVector().getZ());
					}
				}*/
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 * Same scene as testSolve14 but with 1 more particle - runs 254 steps then crashes?
	 * NOTE: Why adding 1 particle to the same scene makes it run longer?
	 */
	@Test
	public void testSolve15()
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
				
				/*if(cycles == 253){
					for(SPHParticle p : ((SPHModel) models.get(0)).getParticles()){
						if(p.getPositionVector().getP() == SPHConstants.LIQUID_TYPE)
						System.out.println(p.getPositionVector().getX() + " / " + p.getPositionVector().getY() + " / " + p.getPositionVector().getZ());
					}
				}*/
			}
		}
		catch (Exception e)
		{
			fail(e.getMessage());
		}

	}
	
	/*
	 * 296 boundary particles + 216 liquid particles (total of 512 particles) with random position and velocity - runs 1 step and then crashes
	 */
	@Test
	public void testSolve216()
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
				
				/*if(cycles == 0){
					for(SPHParticle p : ((SPHModel) models.get(0)).getParticles()){
						if(p.getPositionVector().getP() == SPHConstants.LIQUID_TYPE)
						System.out.println(p.getPositionVector().getX() + " / " + p.getPositionVector().getY() + " / " + p.getPositionVector().getZ());
					}
				}*/
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
	public void testSolvePureLiquidSceneParticlesMove()
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
}
