package org.geppetto.solver.sph;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.runtime.ACompositeNode;
import org.geppetto.core.model.runtime.AspectSubTreeNode;
import org.geppetto.core.model.state.visitors.RemoveTimeStepsVisitor;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.model.sph.x.Vector3DX;
import org.geppetto.solver.sph.internal.FindNaNVisitor;

public class PCISPHTestUtilities {
	public static String readFile(String path) throws IOException
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
	public static void checkStateTreeForNaN(ACompositeNode set, boolean expected)
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
	 * Removes the first state from state tree values
	 * */
	public static void removeFirstStateFromTree(AspectSubTreeNode set)
	{
		RemoveTimeStepsVisitor removeTimeStepVisitor= new RemoveTimeStepsVisitor(1);
		set.apply(removeTimeStepVisitor);
	}
	
	public static void checkModelForOverlappingParticles(SPHModelX model, boolean expected)
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
	 * Counts how many non-boundary particles are in the model.
	 * */
	public static int countNonBoundaryParticles(SPHModelX model)
	{
		int count = 0;
		for(int i = 0; i < model.getNumberOfParticles(); i++)
		{
			Vector3DX positionVector = (Vector3DX) model.getParticles().get(i).getPositionVector();
			
			if (positionVector.getP() != SPHConstants.BOUNDARY_TYPE){
				count++;
			}
		}
		
		return count;
	}
}
