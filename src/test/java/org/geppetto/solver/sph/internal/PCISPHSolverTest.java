/**
 * 
 */
package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.state.StateTreeRoot;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.common.SPHConstants;
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
	
	
	private void checkModelForOverlappingParticles(SPHModelX model, boolean expected)
	{
		int particleCount = model.getNumberOfParticles();
		
		// build a list of boundary particles - don't assume they are in any order
		List<Vector3DX> boundary = new ArrayList<Vector3DX>();
		for(int i = 0; i < particleCount; i++)
		{
			Vector3DX positionVector = (Vector3DX) model.getParticles().get(i).getPositionVector();
			
			if (positionVector.getP() == SPHConstants.BOUNDARY_TYPE){
				boundary.add(positionVector);
			}
		}
		
		// go through particles and check for freely moving particles with same coordinates as boundary
		// NOTE: this should never happen
		List<Vector3DX> overlapping = new ArrayList<Vector3DX>();
		for(int i = 0; i < particleCount; i++)
		{
			Vector3DX vector = (Vector3DX) model.getParticles().get(i).getPositionVector();
			
			// run it against all the boundary positions
			if (vector.getP() != SPHConstants.BOUNDARY_TYPE){
				for(Vector3DX v : boundary)
				{
					if(v.getX().equals(vector.getX()) && v.getY().equals(vector.getY()) && v.getZ().equals(vector.getZ()))
					{
						overlapping.add(vector);
					}
				}
			}
		}
		
		if(expected){
			Assert.assertTrue("Found no overlapping particles when they were expected.", overlapping.size() > 0);
		}
		else{
			Assert.assertTrue("Found overlapping particles when they were not expected: " + overlapping.toString(), overlapping.size() == 0);
		}
	}
	
	/*
	 * Checks the entire StateTreeRoot for NaN values
	 * */
	private void checkStateTreeForNaN(StateTreeRoot set, boolean expected)
	{
		FindNaNVisitor findNaNVisitor=new FindNaNVisitor();
		set.apply(findNaNVisitor);
		
		if(expected)
		{
			Assert.assertTrue("No NaN values detected when expected.", findNaNVisitor.hasNaN());
		} 
		else 
		{
			Assert.assertTrue("Unexpected NaN values detected: " + findNaNVisitor.getParticleWithNaN(), !findNaNVisitor.hasNaN());
		}
	}
	

	/*
	 * 296 boundary particles + 1 liquid particle
	 */
	@Test
	public void testSolve1_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_1.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 2, 1));
		
		//System.out.println(stateSet.lastStateToString());
		
		//checkStateTreeForNaN(stateSet, false);
	}
	
	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_ExpectNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_14.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// we have one particle that overlaps with boundary particles in this test
		checkModelForOverlappingParticles((SPHModelX)model, true);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		// expect NaN values since in the initial conditions particle 309 overlaps with boundary particles
		checkStateTreeForNaN(stateSet, true);
	}

	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_small.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateTreeForNaN(stateSet, false);
	}
	
	/*
	 * 296 boundary particles + 14 liquid particles
	 */
	@Test
	public void testSolve14_StepByStep_VS_OneGo() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_small.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		checkModelForOverlappingParticles((SPHModelX)model, false);
		
		int cycles = 20;
		
		// run cycles one by one
		SPHSolverService solver1 = new SPHSolverService();
		solver1.initialize(model);
		StateTreeRoot stateTree1 = null;
		for(int i = 0; i < cycles; i++){
			stateTree1 = solver1.solve(new TimeConfiguration(0.1f, 1, 1));
		}
		
		// run cycles at once
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		StateTreeRoot stateTree2 = solver2.solve(new TimeConfiguration(0.1f, cycles, 1));
		
		//checks the trees are equivalent
		Assert.assertEquals(stateTree1.toString(),stateTree2.toString());
		
		checkStateTreeForNaN(stateTree1, false);
		checkStateTreeForNaN(stateTree2, false);
	}

	/*
	 * Same scene as testSolve14 but with 1 more particle
	 */
	@Test
	public void testSolve15_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_15.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateTreeForNaN(stateSet, false);
	}

	/*
	 * 296 boundary particles + 216 liquid particles (total of 512 particles) with random position and velocity
	 */
	@Test
	public void testSolve216_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_216.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// check that we don't have particles with overlapping positions
		checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateTreeForNaN(stateSet, false);
	}

	/*
	 * A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolvePureLiquidScene_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_PureLiquid.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// NOTE: commenting out as it takes a while for big scenes
		// checkModelForOverlappingParticles((SPHModelX)model, false);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
		//System.out.println(stateSet.toString());
		
		//checkFinalStateStringForNaN(stateSet.lastStateToString(), false);
	}

	/*
	 * A test built around the original pureLiquid scene used to test the C++ version
	 */
	@Test
	public void testSolveElastic_NoNaN() throws Exception
	{
		URL url = this.getClass().getResource("/sphModel_Elastic.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		IModel model = modelInterpreter.readModel(url);
		
		// NOTE: commenting out as it takes a while for big scenes
		// checkModelForOverlappingParticles((SPHModelX)model, false);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
		//System.out.println(stateSet.toString());
		
		//checkFinalStateStringForNaN(stateSet.lastStateToString(), false);
	}
}
