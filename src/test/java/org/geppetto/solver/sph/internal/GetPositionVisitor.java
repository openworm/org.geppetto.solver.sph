package org.geppetto.solver.sph.internal;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.geppetto.core.model.runtime.VariableNode;
import org.geppetto.core.model.state.visitors.DefaultStateVisitor;
import org.geppetto.model.sph.common.SPHConstants;
import org.geppetto.model.sph.x.Vector3DX;

public class GetPositionVisitor extends DefaultStateVisitor {
	
	private final Integer particleId = 0;
	List<Vector3DX> positionList = new ArrayList<Vector3DX>();
	Vector3DX currentPos = new Vector3DX();
	@Override
	public boolean visitVariableNode(VariableNode node)
	{
		Assert.assertFalse(node.getTimeSeries().size() == 0);
		//positionList.clear();
		int i=0;
		switch(node.getName()){
			case "x":
				currentPos.setX(Float.parseFloat(node.getTimeSeries().get(node.getTimeSeries().size() - 1).toString()));
				break;
			case "y":
				currentPos.setY(Float.parseFloat(node.getTimeSeries().get(node.getTimeSeries().size() - 1).toString()));
				break;
			case "z":
				currentPos.setZ(Float.parseFloat(node.getTimeSeries().get(node.getTimeSeries().size() - 1).toString()));
				break;
			case "p":
				currentPos.setP(Float.parseFloat(node.getTimeSeries().get(node.getTimeSeries().size() - 1).toString()));
				break;
		}
		if(currentPos.getX() != null && currentPos.getY() != null && currentPos.getZ() != null && currentPos.getP() != null){
			if(currentPos.getP() != SPHConstants.BOUNDARY_TYPE){
				positionList.add(currentPos);
			}
			currentPos = new Vector3DX();
		}
		return super.visitVariableNode(node);
	}

	/*public Vector3D getPArticlePosition(int id){
		//return	new Vector3D(0,0,0)	
	}*/
	public Vector3DX getParticlePosition(){
		Vector3DX centerOfMass = new Vector3DX(0.0f,0.0f,0.0f);
		for(Vector3DX v: positionList){
				centerOfMass = Vector3DX.addition(v, centerOfMass);	
		}
		centerOfMass = Vector3DX.mutiplicationOnScalar(centerOfMass, 1.0f/(float)positionList.size());
		positionList.clear();
		return centerOfMass;
	}
	public Vector3DX getParticlePosition(int particleId){
		Vector3DX resultV = positionList.get(particleId);
		positionList.clear();
		return resultV;
	}
}
