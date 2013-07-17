package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.state.StateTreeRoot;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.Vector3D;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.SPHParticleX;
import org.geppetto.solver.sph.KernelsEnum;
import org.geppetto.solver.sph.PCISPHCheckPoint;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class PhysicsLawsTest {

	@Test
	public void testGravitation_SingleParticle() throws Exception {
		// load Java generated scene
		URL url = this.getClass().getResource("/single-elastic-particle.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url);
		
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		
		SPHParticleX p = null;
		for(int i=0; i<model.getNumberOfParticles(); i++){
			if(model.getParticles().get(i).getPositionVector().getP() == SPHConstants.LIQUID_TYPE)
			{
				p = (SPHParticleX)model.getParticles().get(i);
				break;
			}
		}
		
		int steps = 10;
		
		// velocity at time t0
		Vector3D v0 = p.getVelocityVector();
		
		// velocity at time t1
		Vector3D v1 = new Vector3D();
		float t = steps*SPHConstants.TIME_STEP;
		v1.setX((float) (v0.getX() + SPHConstants.GRAVITY_X*t) * SPHConstants.SIMULATION_SCALE);
		v1.setY((float) (v0.getY() + SPHConstants.GRAVITY_Y*t) * SPHConstants.SIMULATION_SCALE);
		v1.setZ((float) (v0.getZ() + SPHConstants.GRAVITY_Z*t) * SPHConstants.SIMULATION_SCALE);
		
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(null, steps, null));
		
		PCISPHCheckPoint check = solver.getCheckpointsMap().get(KernelsEnum.INTEGRATE);
		List<Float> velocity = check.velocity;
		Vector3D v2 = new Vector3D();
		v2.setX(velocity.get(0));
		v2.setY(velocity.get(1));
		v2.setZ(velocity.get(2));
		
		Assert.assertTrue(v1.getX().floatValue() == v2.getX().floatValue());
		Assert.assertTrue(v1.getZ().floatValue() == v2.getZ().floatValue());
		Assert.assertTrue("Y coordinates do not match v1: " + v1.getY().floatValue() + " v2: " + v2.getY().floatValue(), v1.getY().floatValue() == v2.getY().floatValue());
	}

}
