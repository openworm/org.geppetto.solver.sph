/**
 * 
 */
package org.geppetto.solver.sph.internal;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.geppetto.core.model.IModel;
import org.geppetto.core.model.StateInstancePath;
import org.geppetto.core.model.StateSet;
import org.geppetto.core.model.values.AValue;
import org.geppetto.core.model.values.FloatValue;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.SPHModel;
import org.geppetto.model.sph.SPHParticle;
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
	private static final String NAN = "NaN";
	
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
	 * Checks the entire StateSet for NaN values
	 * NOTE: this takes a lot of time if we have many particles - use checkFinalStateStringForNaN for big scenes
	 * */
	private void checkStateSetForNaN(StateSet set, boolean expected)
	{
		List<String> matches = new ArrayList<String>();
		
		Map<StateInstancePath, List<AValue>> stateMap = set.getStatesMap();
		
		for(StateInstancePath k : stateMap.keySet())
		{
			for(AValue v : stateMap.get(k))
			{
				if(v.getStringValue().equals(NAN))
				{
					matches.add(k.toString());
				}
			}
		}
		
		if(expected)
		{
			Assert.assertTrue("No NaN values detected when expected.", matches.size() > 0);
		} 
		else 
		{
			String details = matches.toString();
			Assert.assertTrue("Unexpected NaN values detected: " + details, matches.size() == 0);
		}
	}
	
	/*
	 * Once we get NaN values they will stay like that - so we can just look at the final values
	 * */
	private void checkFinalStateStringForNaN(String stateString, boolean expected)
	{
		if(expected)
		{
			Assert.assertTrue("No NaN values detected when expected.", stateString.contains(NAN));
		} 
		else 
		{
			Assert.assertFalse("Unexpected NaN values detected.", stateString.contains(NAN));
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
		
		SPHSolverService solver = new SPHSolverService(SPHConstants.CPU_PROFILE);
		solver.initialize(model);
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 2, 1));
		
		System.out.println(stateSet.lastStateToString());
		
		checkStateSetForNaN(stateSet, true);
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
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 10, 1));
		
		//System.out.println(stateSet.toString());
		
		// expect NaN values since in the initial conditions particle 309 overlaps with boundary particles
		checkStateSetForNaN(stateSet, true);
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
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateSetForNaN(stateSet, false);
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
		StateSet stateSet1 = null;
		for(int i = 0; i < cycles; i++){
			stateSet1 = solver1.solve(new TimeConfiguration(0.1f, 1, 1));
		}
		
		// run cycles at once
		SPHSolverService solver2 = new SPHSolverService();
		solver2.initialize(model);
		StateSet stateSet2 = solver2.solve(new TimeConfiguration(0.1f, cycles, 1));
		
		// assert state sets equality at the last cycle
		Map<StateInstancePath, List<AValue>> stateMap1 = stateSet1.getStatesMap();
		Map<StateInstancePath, List<AValue>> stateMap2 = stateSet2.getStatesMap();
		
		for(StateInstancePath k : stateMap2.keySet())
		{
			// last cycle of stateMap2
			AValue v2 = stateMap2.get(k).get(cycles - 1);
			// only cycle - the last - of stateMap1
			AValue v1 = stateMap1.get(k).get(0);
			
			Assert.assertTrue(v1.getStringValue().equals(v2.getStringValue()));
		}
		
		checkStateSetForNaN(stateSet1, false);
		checkStateSetForNaN(stateSet2, false);
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
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 20, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateSetForNaN(stateSet, false);
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
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 3, 1));
		
		//System.out.println(stateSet.toString());
		
		checkStateSetForNaN(stateSet, false);
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
		
		// we have one particle that overlaps with boundary particles in this test
		// NOTE: commenting out as it takes a while for big scenes
		// checkModelForOverlappingParticles((SPHModelX)model, true);
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
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
		
		// we have one particle that overlaps with boundary particles in this test
		// NOTE: commenting out as it takes a while for big scenes
		// checkModelForOverlappingParticles((SPHModelX)model, true);

		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		StateSet stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
		
		//System.out.println(stateSet.toString());
		
		//checkFinalStateStringForNaN(stateSet.lastStateToString(), false);
	}
}
