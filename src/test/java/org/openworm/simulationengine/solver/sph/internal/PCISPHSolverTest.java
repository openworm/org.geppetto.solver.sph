/**
 * 
 */
package org.openworm.simulationengine.solver.sph.internal;

import static org.junit.Assert.*;

import java.net.URL;
import java.util.List;

import org.junit.Test;
import org.openworm.simulationengine.core.model.IModel;
import org.openworm.simulationengine.model.sph.SPHModel;
import org.openworm.simulationengine.model.sph.SPHParticle;
import org.openworm.simulationengine.model.sph.common.SPHConstants;
import org.openworm.simulationengine.model.sph.services.SPHModelInterpreterService;
import org.openworm.simulationengine.solver.sph.SPHSolverService;

/**
 * @author matteocantarelli
 * 
 */
public class PCISPHSolverTest
{
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
}
