package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.state.StateTreeRoot;
import org.geppetto.core.model.values.FloatValue;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.Vector3D;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.SPHParticleX;
import org.geppetto.model.sph.x.Vector3DX;
import org.geppetto.solver.sph.KernelsEnum;
import org.geppetto.solver.sph.PCISPHCheckPoint;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

import com.android.dx.cf.code.Simulator;

public class PhysicsLawsTest {

	@Test
	public void testGravitation_SingleParticle() throws Exception {
		// load Java generated scene
		URL url = this.getClass().getResource("/single-elastic-particle.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		solver.setRecordCheckpoint(true);
		
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
		v1.setX((float) (v0.getX() + SPHConstants.GRAVITY_X*t));
		v1.setY((float) (v0.getY() + SPHConstants.GRAVITY_Y*t));
		v1.setZ((float) (v0.getZ() + SPHConstants.GRAVITY_Z*t));
		
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
	
	@Test
	public void test_testGravitation_SingleParticle_1() throws Exception{
		int iterationStart = 0;
		int iterationStop = iterationStart + 100;
		float destination_dis =  0.12f;
		float s = 0.f;
		float s_b = 0.f;
		// load Java generated scene
		URL url = this.getClass().getResource("/gravity_test_one_particle.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		SPHSolverService solver = new SPHSolverService(false);
		solver.initialize(model);
		
		SPHParticleX p = null;
		for(int i=0; i<model.getNumberOfParticles(); i++){
			if(model.getParticles().get(i).getPositionVector().getP() == SPHConstants.ELASTIC_TYPE)
			{
				p = (SPHParticleX)model.getParticles().get(i);
				break;
			}
		}
		Vector3DX initPosition = new Vector3DX(p.getPositionVector());
		int steps = 1;
		Vector3DX end_p;
		GetPositionVisitor positionVisitor = new GetPositionVisitor();
		StateTreeRoot stateSet;
		while(s * SPHConstants.SIMULATION_SCALE <= destination_dis){
			s_b = s;
			stateSet = solver.solve(new TimeConfiguration(null, steps, null));
			stateSet.apply(positionVisitor);
			//stateSet.
			//PCISPHCheckPoint check = solver.getCheckpointsMap().get(KernelsEnum.INTEGRATE);
			//List<Float> position = check.position;
			end_p = positionVisitor.getParticlePosition();
			s = Vector3DX.subtraction(end_p, initPosition).length();
			iterationStart++;
		}
		System.out.println("Distance :" + String.valueOf(s * SPHConstants.SIMULATION_SCALE ) + " m");
		System.out.println("Time :" + String.valueOf( iterationStart * SPHConstants.TIME_STEP) + " s");
	}
	@Test
	public void test_testGravity_Elastic_Cube() throws Exception{
		int iterationStart = 0;
		int iterationStop = iterationStart + 100;
		float destination_dis =  0.02f;
		float s = 0.f;
		float s_b = 0.f;
		URL url = this.getClass().getResource("/gravity_test_many_particle.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		
		SPHSolverService solver = new SPHSolverService(false);
		solver.initialize(model);
		Vector3DX initPosition = new Vector3DX(0.f,0.f,0.f);
		SPHParticleX p = null;
		int particleCount = -1;
		for(int i=0; i<model.getNumberOfParticles(); i++){
			if(model.getParticles().get(i).getPositionVector().getP() == SPHConstants.ELASTIC_TYPE)
			{
				p = (SPHParticleX)model.getParticles().get(i);
				initPosition = Vector3DX.addition(initPosition, new Vector3DX(p.getPositionVector()));
				particleCount++;
			}
		}
		initPosition = Vector3DX.mutiplicationOnScalar(initPosition, 1.f/(float)particleCount);
		int steps = 1;
		Vector3DX end_p;
		GetPositionVisitor positionVisitor = new GetPositionVisitor();
		StateTreeRoot stateSet;
		while(s * SPHConstants.SIMULATION_SCALE <= destination_dis){
			stateSet = solver.solve(new TimeConfiguration(null, steps, null));
			stateSet.apply(positionVisitor);
			end_p = positionVisitor.getParticlePosition();
			s = Vector3DX.subtraction(end_p, initPosition).length();
			iterationStart++;
		}
		System.out.println("Distance :" + String.valueOf(s * SPHConstants.SIMULATION_SCALE ) + " m");
		System.out.println("Time :" + String.valueOf( iterationStart * SPHConstants.TIME_STEP) + " s");
	}
	@Test
	public void test_testFriction_Elastic_Cube() throws Exception{
		int iterationNum = 100;
		int iterationStart = iterationNum - 100;
		int particleId = 100;
		URL url = this.getClass().getResource("/friction_test_1.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url,null,"");
		SPHSolverService solver = new SPHSolverService(false);
		solver.initialize(model);
		Vector3DX initPosition = new Vector3DX(0.f,0.f,0.f);
		int i = 1;
		int steps = 100;
		GetPositionVisitor positionVisitor = new GetPositionVisitor();
		StateTreeRoot stateSet;
		SPHParticleX p = null;
		int particleCount = 0;
		for(int j=0; j<model.getNumberOfParticles(); j++){
			if(model.getParticles().get(j).getPositionVector().getP() == SPHConstants.ELASTIC_TYPE)
			{
				p = (SPHParticleX)model.getParticles().get(j);
				initPosition = Vector3DX.addition(initPosition, new Vector3DX(p.getPositionVector()));
				particleCount++;
			}
		}
		initPosition = Vector3DX.mutiplicationOnScalar(initPosition, 1.f/(float)particleCount);
		Vector3DX end_p;
		stateSet = solver.solve(new TimeConfiguration(null, steps, null));
		stateSet.apply(positionVisitor);
		end_p = positionVisitor.getParticlePosition();
		float result = Vector3DX.subtraction(initPosition, end_p).length();
		System.out.println("Distance :" + String.valueOf(result * SPHConstants.SIMULATION_SCALE ) + " m");
		System.out.println("Time :" + String.valueOf( (iterationNum - iterationStart) * SPHConstants.TIME_STEP) + " s");
	}
}

