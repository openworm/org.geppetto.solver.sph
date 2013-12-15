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
package org.geppetto.solver.sph;

import java.util.List;

import org.bridj.Pointer;
import org.geppetto.core.model.state.CompositeStateNode;
import org.geppetto.core.model.state.SimpleStateNode;
import org.geppetto.core.model.state.visitors.DefaultStateVisitor;
import org.geppetto.core.model.values.FloatValue;
import org.geppetto.core.model.values.ValuesFactory;

public class UpdateSPHWatchTreeVisitor extends DefaultStateVisitor {

	private FloatValue _xP, _yP, _zP, _pP;
	private FloatValue _xV, _yV, _zV, _pV;
	private Pointer<Float> _positionPtr;
	private Pointer<Float> _velocityPtr;
	private List<String> _watchlist;

	public UpdateSPHWatchTreeVisitor(Pointer<Float> positionPtr, Pointer<Float> velocityPtr, List<String> watchList)
	{
		_positionPtr = positionPtr;
		_velocityPtr = velocityPtr;
		_watchlist = watchList;
	}

	@Override
	public boolean inCompositeStateNode(CompositeStateNode node)
	{
		if(node.isArray())
		{
			int index = node.getIndex()*4;

			_xP = ValuesFactory.getFloatValue(_positionPtr.get(index));
			_yP = ValuesFactory.getFloatValue(_positionPtr.get(index + 1));
			_zP = ValuesFactory.getFloatValue(_positionPtr.get(index + 2));
			_pP = ValuesFactory.getFloatValue(_positionPtr.get(index + 3));
			
			_xV = ValuesFactory.getFloatValue(_velocityPtr.get(index));
			_yV = ValuesFactory.getFloatValue(_velocityPtr.get(index + 1));
			_zV = ValuesFactory.getFloatValue(_velocityPtr.get(index + 2));
			_pV = ValuesFactory.getFloatValue(_velocityPtr.get(index + 3));
				
		}
		return super.inCompositeStateNode(node);
	}

	@Override
	public boolean outCompositeStateNode(CompositeStateNode node)
	{
		if(node.isArray())
		{
			_xV = _yV = _zV = _pV = null;
			_xP = _yP = _zP = _pP = null;
		}
		return super.outCompositeStateNode(node);
	}

	@Override
	public boolean visitSimpleStateNode(SimpleStateNode node)
	{
		// TODO: 1. check if this is the node we are watching
		// TODO: 2. get value of interest
		// TODO: 3. node.addValue
		
		return super.visitSimpleStateNode(node);
	}

}
