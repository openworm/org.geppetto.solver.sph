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

public enum KernelsEnum {
	    CLEAR_BUFFERS("clearBuffers"),
	    HASH_PARTICLES("hashParticles"),
	    SORT("sort"),
	    SORT_POST_PASS("sortPostPass"),
	    INDEX("indexx"),
	    INDEX_POST_PASS("indexxPostPass"),
	    FIND_NEIGHBORS("findNeighbors"),
	    COMPUTE_DENSITY("pcisph_computeDensity"),
	    COMPUTE_FORCES_INIT_PRESSURE("pcisph_computeForcesAndInitPressure"),
	    COMPUTE_ELASTIC_FORCES("pcisph_computeElasticForces"),
	    PREDICT_POSITION("pcisph_predictPositions"),
	    PREDICT_DENSITY("pcisph_predictDensity"),
	    CORRECT_PRESSURE("pcisph_correctPressure"),
	    COMPUTE_PRESSURE_FORCE_ACCELERATION("pcisph_computePressureForceAcceleration"),
	    PREDICTIVE_LOOP("predictiveLoop"),
	    INTEGRATE("pcisph_integrate"),
	    ;
	    
	    private KernelsEnum(final String text) {
	        this.text = text;
	    }

	    private final String text;

	    @Override
	    public String toString() {
	        return text;
	    }
}
