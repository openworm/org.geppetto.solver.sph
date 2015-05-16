/*******************************************************************************
 * The MIT License (MIT)
 * 
 * Copyright (c) 2011 - 2015 OpenWorm.
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

import java.math.BigDecimal;
import java.util.HashSet;
import java.util.Set;

import org.geppetto.core.model.quantities.PhysicalQuantity;
import org.geppetto.core.model.runtime.CompositeNode;
import org.geppetto.core.model.runtime.VariableNode;
import org.geppetto.core.model.state.visitors.DefaultStateVisitor;
import org.geppetto.core.model.values.AValue;
import org.geppetto.model.sph.Vector3D;

/**
 * @author giovanni@openworm.org
 */
public class CompareStateVisitor extends DefaultStateVisitor
{
	private Integer currentID = null;
	private String[] referenceStates = null;
	private Set<Integer> mismatchingIDs = new HashSet<>();
	
	private final String X = "x";
	private final String Y = "y";
	private final String Z = "z";
	private final String P = "p";
	
	public CompareStateVisitor(String[] referenceStates) {
		super();
		this.referenceStates = referenceStates;
	}

	public Set<Integer> getMismatches()
	{
		return mismatchingIDs;
	}
	
	@Override
	public boolean inCompositeNode(CompositeNode node) {
		if(node.isArray())
			currentID = node.getIndex();
		
		return super.inCompositeNode(node);
	}
	
	@Override
	public boolean outCompositeNode(CompositeNode node) {
		if(node.isArray())
			currentID = null;
		
		return super.inCompositeNode(node);
	}

	@Override
	public boolean visitVariableNode(VariableNode node)
	{
		// get last step
		PhysicalQuantity p = node.getTimeSeries().get(node.getTimeSeries().size() - 1);
		AValue v = p.getValue();
		Float nodeVal = Float.parseFloat(v.getStringValue());
		
		String refValues = referenceStates[currentID];
		Vector3D refVector = get3DVector(refValues);
		Float refVal = null;
		
		String name = node.getName();
		switch (name) {
        case X:  
        	refVal = refVector.getX();
            break;
        case Y:
        	refVal = refVector.getY();
        	break;
        case Z:
        	refVal = refVector.getZ();
        	break;
        case P:
        	refVal = refVector.getP();
        	break;
		}
		
		// round to 3rd decimal digit for X - Y - Z
		if((name.equals(X) || name.equals(Y) || name.equals(Z)) &&
		   !(round(nodeVal.floatValue(), 0) == round(refVal, 0)))
		{
			// add ID to mismatching set
			mismatchingIDs.add(currentID);
		}
		else if (name.equals(P) && 
				 !(round(nodeVal.floatValue(), 0) == round(refVal, 0)))
		{
			// add ID to mismatching set
			mismatchingIDs.add(currentID);
		}
		
		return super.visitVariableNode(node);
	}
	
	private float round(float d, int decimalPlace) 
	{
		BigDecimal bd = new BigDecimal(Float.toString(d));
		bd = bd.setScale(decimalPlace, BigDecimal.ROUND_HALF_UP);
		return bd.floatValue();
	}
	
	private Vector3D get3DVector(String triplet)
	{
		Vector3D v = new Vector3D();
		String[] coordinates = triplet.split("\t");
		v.setX(new Float(coordinates[0].trim()));
		v.setY(new Float(coordinates[1].trim()));
		v.setZ(new Float(coordinates[2].trim()));
		if (coordinates.length > 3)
		{
			v.setP(new Float(coordinates[3].trim()));
		}
		return v;
	}
}
