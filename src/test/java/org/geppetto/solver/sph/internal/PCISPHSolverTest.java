/**
 * 
 */
package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.StateSet;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.SPHModel;
import org.geppetto.model.sph.SPHParticle;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.Vector3DX;
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
					// non boundary particles for these tests are never on the origin after one cycle if everything works
					Assert.assertTrue("Something is not working, position of the first non boundary particle shouldn't be the origin. Cycles executed:" + (cycles + 1), !(checkX && checkY && checkZ));
					break;
				}
			}
		}
	}

	private void printCoordinates(List<IModel> models)
	{
		SPHModelX mod = (SPHModelX) models.get(0);

		int i = 1;
		for(SPHParticle p : mod.getParticles())
		{
			Vector3DX pos = (Vector3DX) p.getPositionVector();
			Vector3DX vel = (Vector3DX) p.getVelocityVector();

			System.out.println("#" + i + " position x:" + pos.getX() + " y:" + pos.getY() + " z:" + pos.getZ() + " p:" + pos.getP());
			System.out.println("#" + i + " velocity x:" + vel.getX() + " y:" + vel.getY() + " z:" + vel.getZ() + " p:" + vel.getP());
			i++;
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
		URL url = new URL("https://www.dropbox.com/s/eshuozw196k3vci/sphModel_14.xml?dl=1");

		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 10, 1));
		System.out.println(stateSet.toString());
	}

	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoCrash2() throws Exception
	{
		SPHSolverService solver = new SPHSolverService();
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		URL url = new URL("https://www.dropbox.com/s/8869zlz971ogyra/sphModel_small.xml?dl=1");

		// TODO: values go to NaN around 19 steps - figure out why!
		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 19, 1));
		System.out.println(stateSet.toString());
	}

	/*
	 * Same scene as testSolve14 but with 1 more particle
	 */
	@Test
	public void testSolve15_NoCrash() throws Exception
	{
		SPHSolverService solver = new SPHSolverService();
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		URL url = new URL("https://www.dropbox.com/s/9kx2p8qspdgphd4/sphModel_15.xml?dl=1");
		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 10, 1));
		System.out.println(stateSet.toString());
	}

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
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 10, 1));
		System.out.println(stateSet.toString());
	}

	/*
	 * A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolvePureLiquidScene_ParticlesMoving() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");

		SPHSolverService solver = new SPHSolverService();
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();

		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 10, 1));
		System.out.println(stateSet.toString());
	}

	/*
	 * A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolveElastic_NoCrash() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");

		SPHSolverService solver = new SPHSolverService();
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();

		IModel model = modelInterpreter.readModel(url);
		solver.initialize(model);
		StateSet stateSet=solver.solve(new TimeConfiguration(0.1f, 10, 1));
		System.out.println(stateSet.toString());
	}
}
