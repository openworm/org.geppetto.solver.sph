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
import org.geppetto.core.model.typesystem.values.ColladaValue;
import org.geppetto.core.model.typesystem.values.CompositeValue;
import org.geppetto.core.model.typesystem.values.Cylinder;
import org.geppetto.core.model.typesystem.values.FloatValue;
import org.geppetto.core.model.typesystem.values.ParticleValue;
import org.geppetto.core.model.typesystem.values.SphereValue;
import org.geppetto.core.model.typesystem.values.ValuesFactory;
import org.geppetto.core.model.typesystem.visitor.AnalysisVisitor;
import org.geppetto.core.visualisation.model.Point;

/**
 * @author matteocantarelli
 * 
 * This method updates the particles already present in the tree
 * adding new values as found on the position pointer
 */
public class UpdateSPHVisualizationTreeVisitor extends AnalysisVisitor
{

	private FloatValue _xV, _yV, _zV, _pV;
	private Pointer<Float> _positionPtr;
	
	public UpdateSPHVisualizationTreeVisitor(Pointer<Float> positionPtr)
	{
		this._positionPtr = positionPtr;
	}

	@Override
	public boolean inCompositeNode(CompositeValue node)
	{		
		return super.inCompositeNode(node);
	}
	
	@Override
	public boolean outCompositeNode(CompositeValue node)
	{
		return super.outCompositeNode(node);
	}
	
	@Override
	public boolean visitSphereNode(SphereValue node){

		if(node.getPosition() != null){
			Point newPosition = new Point();
			double x = this._xV.getAsDouble();
			double y = this._yV.getAsDouble();
			double z = this._zV.getAsDouble();
			
			if(!Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)){
				newPosition.setX(x);
				newPosition.setY(y);
				newPosition.setZ(z);
				node.setPosition(newPosition);
			}
		}
		return super.visitSphereNode(node);
	}
	
	@Override
	public boolean visitCylinderNode(Cylinder node){

		if(node.getPosition() != null){
			Point newPosition = new Point();
			double x = this._xV.getAsDouble();
			double y = this._yV.getAsDouble();
			double z = this._zV.getAsDouble();
			
			if(!Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)){
				newPosition.setX(x);
				newPosition.setY(y);
				newPosition.setZ(z);
				node.setPosition(newPosition);
			}
		}
		return super.visitCylinderNode(node);
	}
	
	@Override
	public boolean visitColladaNode(ColladaValue node){
		if(node.getPosition() != null){
			Point newPosition = new Point();
			double x = this._xV.getAsDouble();
			double y = this._yV.getAsDouble();
			double z = this._zV.getAsDouble();
			
			if(!Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)){
				newPosition.setX(x);
				newPosition.setY(y);
				newPosition.setZ(z);
				node.setPosition(newPosition);
			}
		}
		return super.visitColladaNode(node);
	}
	
	@Override
	public boolean visitParticleNode(ParticleValue node){
		if(node.getPosition() != null){
			
			int index = node.getIndex() * 4;
			_xV = ValuesFactory.getFloatValue(_positionPtr.get(index));
			_yV = ValuesFactory.getFloatValue(_positionPtr.get(index + 1));
			_zV = ValuesFactory.getFloatValue(_positionPtr.get(index + 2));
 			_pV = ValuesFactory.getFloatValue(_positionPtr.get(index + 3));
			
 			Point newPosition = new Point();
			double x = this._xV.getAsDouble();
			double y = this._yV.getAsDouble();
			double z = this._zV.getAsDouble();
			float p = this._pV.getAsFloat();
			
			if(!Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)){
				newPosition.setX(x);
				newPosition.setY(y);
				newPosition.setZ(z);
				node.setPosition(newPosition);
			}
			
			if(!Float.isNaN(p)){
				node.setParticleKind(p);
			}
		}
		return super.visitParticleNode(node);
	}
}
