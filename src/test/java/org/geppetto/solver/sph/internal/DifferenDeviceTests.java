/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

package org.geppetto.solver.sph.internal;

import org.junit.Before;
import org.junit.Test;

import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.geppetto.core.model.runtime.AspectNode;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.Membrane;
import org.geppetto.model.sph.Vector3D;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.Vector3DX;
import org.geppetto.solver.sph.BuffersEnum;
import org.geppetto.solver.sph.KernelsEnum;
import org.geppetto.solver.sph.PCISPHCheckPoint;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;

public class DifferenDeviceTests {
	
	AspectNode _sphAspect = null;

	@Before
	public void runBeforeEveryTest()
	{
		_sphAspect=new AspectNode();
	}
	
	@Test
	public void test_different_machine_result() throws Exception{
		/* I plane it will be work in two modes:
		 * 1) First it logging information about data buffers into the files
		 * 2) Second it compares information from logs with info from current run 
		 */
		boolean logInfo = true;
		int iterationCount = 10;
		logOrCompare(logInfo, this.getClass().getResource("/cube_with_membranes_water_inside.xml"), KernelsEnum.INTEGRATE, iterationCount);
		
	}
	/**Logging information or compare simulation results with results 
	 * took from another source. We consider only position for now.
	 * @param logInfo
	 * indicate should test compare or logging data
	 * @param modelURL
	 * @param kernel
	 * which kernel should be logged you can use any kernel from enumeration KernelsEnum if
	 * @param iterationCount
	 * number of iteration 
	 * @throws Exception
	 */
	private void logOrCompare(boolean logInfo, URL modelURL, KernelsEnum checkpoint, int iterationCount) throws Exception
	{
		// load reference values at various steps from logs version
		String path = "/home/serg/git/openworm/org.geppetto.solver.sph/src/test/resources/results/DifferentMachinesTest/";
		Map<BuffersEnum, String[]> checkpointReferenceValuesMap = new HashMap<BuffersEnum, String[]>();
		
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		
		// load Java generated scene
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(modelURL,null,"");
		
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		if(checkpoint !=null ){
			solver.setCheckKernel(checkpoint);
			solver.setCheckBuffer(BuffersEnum.POSITION);
		}
		Map<BuffersEnum, Integer> mismatchingValuesPerBuffers = new LinkedHashMap<BuffersEnum, Integer>();
		// calculate step
		int iteration = 0;
		// get buffer sizes
		Map<BuffersEnum, Integer> dimensions = solver.getBuffersSizeMap();
		Map<BuffersEnum,List<Map<String,Vector3D>>> partcileEvolution = new HashMap<BuffersEnum,List<Map<String,Vector3D>>>();
		while(iteration <= iterationCount){
			if(!logInfo){
				logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/DifferentMachinesTest/position_" + iteration + ".txt"));
				for (Map.Entry<BuffersEnum, URL> entry : logs.entrySet())
				{
					String fileContent = PCISPHTestUtilities.readFile(logs.get(entry.getKey()).getPath());
					checkpointReferenceValuesMap.put(entry.getKey(), fileContent.split(System.getProperty("line.separator")));
				}
			}
			solver.solve(new TimeConfiguration(null, 1, null),_sphAspect);
			// get checkpoint of interest
			PCISPHCheckPoint checkpoint_values = solver.getCheckpointsMap().get(checkpoint);
	      	int j=0;
	      	// compare buffer values at first checkpoint
	      	if(!logInfo){
				for (Map.Entry<BuffersEnum, String[]> entry : checkpointReferenceValuesMap.entrySet())
				{
					switch (entry.getKey()) {
				        case POSITION:
			        		List<Float> pos_calculatedValues = checkpoint_values.position;
					      	int pos_mismatches = 0;
					      	boolean isDifferent = false;
					      	//In file storred information only for moving particles boundary particles isn't taking into account
					      	Assert.assertTrue(solver._numOfElasticP + solver._numOfLiquidP == checkpointReferenceValuesMap.get(BuffersEnum.POSITION).length);
					      	j=0;
					      	for(int i = 0; i < (solver._numOfElasticP + solver._numOfLiquidP)*4; i = i + 4)
					      	{
					      		// get vector from ref values
					      		Vector3D pos_vector = get3DVector(checkpointReferenceValuesMap.get(BuffersEnum.POSITION)[j]);
					      		
					      		// it sucks, but all the if statements are separate to facilitate debugging
					      		if(pos_vector.getX().floatValue() != pos_calculatedValues.get(i).floatValue()) 
					      		{
					      			pos_mismatches++;
					      			isDifferent = true;
					      		}
					      		if(pos_vector.getY().floatValue() != pos_calculatedValues.get(i + 1).floatValue()) 
					      		{
					      			pos_mismatches++;
					      			isDifferent = true;
					      		}
					      		if(pos_vector.getZ().floatValue() != pos_calculatedValues.get(i + 2).floatValue()) 
					      		{
					      			pos_mismatches++;
					      			isDifferent = true;
					      		}
					      		if(pos_vector.getP().floatValue() != pos_calculatedValues.get(i + 3).floatValue()) 
					      		{
					      			pos_mismatches++;
					      			isDifferent = true;
					      		}
					      		if(isDifferent){
					      			Map<String,Vector3D> temp = new HashMap<String, Vector3D>();
					      			temp.put("GEPPETTO_C", new Vector3DX(pos_calculatedValues.get(i).floatValue(),pos_calculatedValues.get(i + 1).floatValue(),pos_calculatedValues.get(i + 2).floatValue()));
					      			temp.put("GEPPETTO_O", pos_vector);
					      			if(partcileEvolution.containsKey(BuffersEnum.POSITION))
					      				partcileEvolution.get(BuffersEnum.POSITION).add(temp);
					      			else{
					      				List<Map<String,Vector3D>> t_list= new ArrayList<Map<String,Vector3D>>();
					      				t_list.add(temp);
					      				partcileEvolution.put(BuffersEnum.POSITION, t_list);
					      			}
					      			isDifferent = false;
					      		}
					      		j++;
					      	}
					      	break;
					    default:
					    	break;
					}
				}
			}
			else{
	      		List<Float> pos_calculatedValues = checkpoint_values.position;
	      		// As I understood I cannot write on a Resource folder (I cannot create file there)
	      		// So I put here absolute path to folder it could be a relative.
	      		PrintWriter writer = new PrintWriter(path + "position_" + iteration + ".txt", "UTF-8");
	      		for(int i = 0; i < dimensions.get(BuffersEnum.POSITION)/2; i = i + 4)
		      	{
	      			if(pos_calculatedValues.get(i + 3) != SPHConstants.BOUNDARY_TYPE){
		      			String s = pos_calculatedValues.get(i) + "\t" + pos_calculatedValues.get(i + 1) + "\t" + pos_calculatedValues.get(i + 2) + "\t" + pos_calculatedValues.get(i + 3);
		      			writer.println(s);
	      			}
	      		}
	      		writer.close();
	      	}
	      	iteration++;
		}
		if(!logInfo){
			int j=0;
			if(partcileEvolution.size() > 0)
				System.out.println("Smulation shows different results comparing with what was found in resources/results/DifferentMachinesTest/ folder");
			else
				System.out.println("Smulation shows the same results comparing with what was found in resources/results/DifferentMachinesTest/ folder");
			for(Map.Entry<BuffersEnum,List<Map<String,Vector3D>>> l:partcileEvolution.entrySet()){
				System.out.println("Diffrents on step " + j);
				for(Map<String,Vector3D> m:l.getValue()){
					System.out.println("GEPPETTO LOGGED X = " + m.get("GEPPETTO_O").getX().floatValue() + " Y = " + m.get("GEPPETTO_O").getY().floatValue() + " Z = " + m.get("GEPPETTO_O").getZ().floatValue());
		  			System.out.println("GEPPETTO CURRENT X = " + m.get("GEPPETTO_C").getX().floatValue() + " Y = " + m.get("GEPPETTO").getY().floatValue() + " Z = " + m.get("GEPPETTO").getZ().floatValue());
				}
				j++;
			}
		}else{
			System.out.println("Information about evolution of position is stored into the files in folder" + path);
		}
	}
	
	private Vector3D get3DVector(String values)
	{
		Vector3D v = new Vector3D();
		String[] coordinates = values.split("\t");
		
		if (coordinates.length > 0) v.setX(new Float(coordinates[0].trim()));
		if (coordinates.length > 1) v.setY(new Float(coordinates[1].trim()));
		if (coordinates.length > 2) v.setZ(new Float(coordinates[2].trim()));
		if (coordinates.length > 3)	v.setP(new Float(coordinates[3].trim()));

		return v;
	}
}
