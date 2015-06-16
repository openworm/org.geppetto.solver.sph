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

import junit.framework.Assert;

import org.geppetto.core.model.quantities.PhysicalQuantity;
import org.geppetto.core.model.quantities.Quantity;
import org.geppetto.core.model.runtime.VariableNode;
import org.geppetto.core.model.state.visitors.DefaultStateVisitor;
import org.geppetto.core.model.values.AValue;

/**
 * @author matteo@openworm.org
 * @author giovanni@openworm.org
 */
public class FindNaNVisitor extends DefaultStateVisitor
{
	private static final String NAN = "NaN";
	
	private boolean _hasNaN=false;
	private String _particleWithNaN=null;
	
	public boolean hasNaN()
	{
		return _hasNaN;
	}

	@Override
	public boolean visitVariableNode(VariableNode node)
	{
		Assert.assertFalse(node.getTimeSeries().size() == 0);
		
		int i=0;
		for(Quantity p:node.getTimeSeries()) 
		{
			AValue v = p.getValue();
			if(v.getStringValue().equals(NAN))	
			{
				doStopVisiting();
				_particleWithNaN=node.getParent().getParent().getName()+"{"+i+"}";
				_hasNaN=true;
				break;
			}
			i++;
		}
		return super.visitVariableNode(node);
	}

	public String getParticleWithNaN()
	{
		return _particleWithNaN;
	}

}
