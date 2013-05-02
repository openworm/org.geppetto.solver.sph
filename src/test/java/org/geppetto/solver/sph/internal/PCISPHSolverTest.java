/**
 * 
 */
package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.StateSet;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.SPHModel;
import org.geppetto.model.sph.SPHParticle;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

/**
 * @author matteo@openworm.org
 * @author giovanni@openworm.org
 */
public class PCISPHSolverTest
{

	private void checkModels(List<IModel> models, int cycles)
	{
		for(IModel m : models)
		{
			SPHModel sphm = (SPHModel) m;
			Assert.assertNotNull(sphm.getParticles());
			Assert.assertTrue(!sphm.getParticles().isEmpty());
			for(SPHParticle p : sphm.getParticles())
			{
				if(!p.getPositionVector().getP().equals(3.1f))
				{
					boolean checkX = p.getPositionVector().getX().equals(0f);
					boolean checkY = p.getPositionVector().getY().equals(0f);
					boolean checkZ = p.getPositionVector().getZ().equals(0f);
					//non boundary particles for these tests are never on the origin after one cycle if everything works
					Assert.assertTrue("Something is not working, position of the first non boundary particle shouldn't be the origin. Cycles executed:"+(cycles+1), !(checkX && checkY && checkZ));
					break;
				}
			}
		}
	}

	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoCrash1() throws Exception
	{
		SPHSolverService solver = new SPHSolverService();
		
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		//URL url = new URL("https://www.dropbox.com/s/eshuozw196k3vci/sphModel_14.xml?dl=1");
		URL url = new URL("https://dl.dropboxusercontent.com/u/7538688/sphModel_14.xml?dl=1");
		
		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 2, 1));
		System.out.println(stateSet.toString());
	}


//	/*
//	 * 296 boundary particles + 14 liquid particles
//	 */
//	@Test
//	public void testSolve14_NoCrash2() throws Exception
//	{
//		SPHSolverService solver = new SPHSolverService();
//		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
//		URL url = new URL("https://www.dropbox.com/s/8869zlz971ogyra/sphModel_small.xml?dl=1");
//		List<IModel> models = modelInterpreter.readModel(url);
//		for(int cycles = 0; cycles < 18; cycles++)
//		{
//			models = solver.solve(models, null).get(0);
//			checkModels(models,cycles);
//		}
//	}
//
//	/*
//	 * Same scene as testSolve14 but with 1 more particle
//	 */
//	@Test
//	public void testSolve15_NoCrash() throws Exception
//	{
//		SPHSolverService solver = new SPHSolverService();
//		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
//		URL url = new URL("https://www.dropbox.com/s/9kx2p8qspdgphd4/sphModel_15.xml?dl=1");
//		List<IModel> models = modelInterpreter.readModel(url);
//		for(int cycles = 0; cycles < 254; cycles++)
//		{
//			models = solver.solve(models, null).get(0);
//			checkModels(models,cycles);
//		}
//	}
//
	/*
	 * 296 boundary particles + 216 liquid particles (total of 512 particles) with random position and velocity
	 */
	@Test
	public void testSolve216_NoCrash() throws Exception
	{
		SPHSolverService solver = new SPHSolverService();
		
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		URL url = new URL("https://www.dropbox.com/s/lerz4rkx75nq0bk/sphModel_216.xml?dl=1");
		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 9, 1));
		System.out.println(stateSet.toString());
		
	}
//
//	/*
//	 * A test built around the original pureLiquid scene used to test the C++ version
//	 */
//	@Test
//	public void testSolvePureLiquidScene_ParticlesMoving() throws Exception
//	{
//		URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");
//
//		SPHSolverService solver = new SPHSolverService();
//		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
//
//		List<IModel> initial_models = modelInterpreter.readModel(url);
//		List<IModel> models = new ArrayList<IModel>(initial_models);
//
//		for(int cycles = 0; cycles < 1; cycles++)
//		{
//			models = solver.solve(models, null).get(0);
//			int pd = ((SPHModelX) initial_models.get(0)).compareTo((SPHModelX) (models.get(0)));
//			System.out.println("Particles different at cycle " + cycles + ": " + pd);
//			Assert.assertFalse(pd == 0);
//			initial_models = models;
//			checkModels(models,cycles);
//		}
//	}
//
//	/*
//	 * A test built around the original pureLiquid scene used to test the C++ version
//	 */
//	@Test
//	public void testSolveElastic_NoCrash() throws Exception
//	{
//		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
//
//		SPHSolverService solver = new SPHSolverService();
//		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
//
//		List<IModel> initial_models = modelInterpreter.readModel(url);
//		List<IModel> models = new ArrayList<IModel>(initial_models);
//
//		for(int cycles = 0; cycles < 2; cycles++)
//		{
//			models = solver.solve(models, null).get(0);
//			checkModels(models, cycles);
//		}
//	}
}
