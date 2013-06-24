package org.geppetto.solver.sph.internal;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import junit.framework.Assert;

import org.geppetto.core.model.state.StateTreeRoot;
import org.geppetto.core.simulation.TimeConfiguration;

import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class StepValidationTest {

	@Test
	public void ValidateLiquidScene_16974() throws Exception {
		// load reference values at various steps from C++ version
		String position1 = readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_1.txt").getPath());
		String position5 = readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_5.txt").getPath());
		String position10 = readFile(StepValidationTest.class.getResource("/results/liquid_16974/position_log_10.txt").getPath());
		
		Map<Integer, String[]> referenceValuesMap = new HashMap<Integer, String[]>();
		referenceValuesMap.put(1, position1.split(System.getProperty("line.separator")));
		referenceValuesMap.put(5, position5.split(System.getProperty("line.separator")));
		referenceValuesMap.put(10, position10.split(System.getProperty("line.separator")));
		
		// load Java generated scene
		URL url = this.getClass().getResource("/sphModel_liquid_16974.xml");
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(url);
		
		// 3. assert number of particles is fine
		for (Map.Entry<Integer, String[]> entry : referenceValuesMap.entrySet())
		{
		    Assert.assertTrue("number of lines on positions and number of particles on sphModel do not match", model.getParticles().size() == entry.getValue().length);
		}
		
		SPHSolverService solver = new SPHSolverService();
		solver.initialize(model);
		
		Map<Integer,Set<Integer>> mismatchingSetsMap = new HashMap<Integer, Set<Integer>>();
		for(int i = 0; i < 10; i++)
		{
			// calculate step
			StateTreeRoot stateSet = solver.solve(new TimeConfiguration(0.1f, 1, 1));
			
			// get reference values
			int step = i + 1;
			String[] referenceValues = referenceValuesMap.get(step);
			
			checkStateTreeForNaN(stateSet, false);
			
			if( referenceValuesMap.containsKey(step) )
			{
				// compare state tree with logged values for each recorded step
				CompareStateVisitor compareVisitor = new CompareStateVisitor(referenceValues);
				stateSet.apply(compareVisitor);
				mismatchingSetsMap.put(step, compareVisitor.getMismatches());
			}
		}
		
		// check mismatching stuff
		for (Map.Entry<Integer,Set<Integer>> entry : mismatchingSetsMap.entrySet())
		{
			Integer step = entry.getKey();
			// assert there are no mismatching values
			Assert.assertTrue(entry.getValue().size() + " mismatching values found on step " + step, entry.getValue().size() == 0);
		}
	}
	
	@Test
	public void ValidateElasticScene_SingleBundle() throws Exception {
		fail("Not yet implemented");
	}
	
	private String readFile(String path) throws IOException
	{
		FileInputStream stream = new FileInputStream(new File(path));
		try
		{
			FileChannel fc = stream.getChannel();
			MappedByteBuffer bb = fc.map(FileChannel.MapMode.READ_ONLY, 0, fc.size());
			/* Instead of using default, pass in a decoder. */
			return Charset.defaultCharset().decode(bb).toString();
		}
		finally
		{
			stream.close();
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
}
