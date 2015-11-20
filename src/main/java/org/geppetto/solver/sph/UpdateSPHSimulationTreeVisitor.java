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
package org.geppetto.solver.sph;

import org.bridj.Pointer;
import org.geppetto.core.model.typesystem.values.FloatValue;
import org.geppetto.core.model.typesystem.values.PhysicalQuantityValue;
import org.geppetto.core.model.typesystem.values.QuantityValue;
import org.geppetto.core.model.typesystem.values.ValuesFactory;
import org.geppetto.core.model.typesystem.values.VariableValue;
import org.geppetto.core.model.typesystem.visitor.AnalysisVisitor;

public class UpdateSPHSimulationTreeVisitor extends AnalysisVisitor {

	private Pointer<Float> _positionPtr;

	public UpdateSPHSimulationTreeVisitor(Pointer<Float> positionPtr)
	{
		_positionPtr = positionPtr;
	}

	@Override
	public boolean visitVariableNode(VariableValue node)
	{
		// 1. figure out which of the variables being watched this node represents
		String fullName = node.getInstancePath();
		
		// 2. get value of interest (Nth particle from relevant results arrays)
		// extract index from string
		Integer particleIndex = null;
		if (fullName.indexOf("[")!=-1) {
			String particleID = fullName.substring(fullName.indexOf("[")+1, fullName.indexOf("]"));
			particleIndex = Integer.parseInt(particleID)*4;
		}
		
		// use index to retrieve values
		FloatValue _xV = ValuesFactory.getFloatValue(_positionPtr.get(particleIndex));
		FloatValue _yV = ValuesFactory.getFloatValue(_positionPtr.get(particleIndex + 1));
		FloatValue _zV = ValuesFactory.getFloatValue(_positionPtr.get(particleIndex + 2));
		
		// 3. node.addValue
		if(node.getId().equals("x"))
		{
			QuantityValue q = new QuantityValue();
			q.setValue(_xV);
			node.addQuantity(q);
		}
		else if(node.getId().equals("y"))
		{
			QuantityValue q = new QuantityValue();
			q.setValue(_yV);
			node.addQuantity(q);
		}
		else if(node.getId().equals("z"))
		{
			QuantityValue q = new QuantityValue();
			q.setValue(_zV);
			node.addQuantity(q);
		}
		
		return super.visitVariableNode(node);
	}

}
