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

import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.geppetto.core.model.runtime.AspectNode;
import org.geppetto.core.simulation.TimeConfiguration;
import org.geppetto.model.sph.Vector3D;
import org.geppetto.model.sph.services.SPHModelInterpreterService;
import org.geppetto.model.sph.x.SPHModelX;
import org.geppetto.solver.sph.BuffersEnum;
import org.geppetto.solver.sph.KernelsEnum;
import org.geppetto.solver.sph.PCISPHCheckPoint;
import org.geppetto.solver.sph.PCISPHTestUtilities;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;

public class StepValidationWithCheckpointsTest {

	@Test
	public void testCheckpoints_780_CLEARBUFFERS() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_density_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_gridcellindex_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_gridcellindexfixedup_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_index_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_indexback_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_neighbormap_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_position_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_pressure_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_sortedposition_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_sortedvelocity_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/01_velocity_log_runClearBuffers_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.CLEAR_BUFFERS, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_HASHPARTICLES() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_density_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_gridcellindex_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_gridcellindexfixedup_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_index_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_indexback_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_neighbormap_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_position_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_pressure_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_sortedposition_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_sortedvelocity_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/02_velocity_log_runHashParticles_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.HASH_PARTICLES, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_SORT() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_density_log_runSort_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_gridcellindex_log_runSort_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_gridcellindexfixedup_log_runSort_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_index_log_runSort_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_indexback_log_runSort_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_neighbormap_log_runSort_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_position_log_runSort_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_pressure_log_runSort_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_sortedposition_log_runSort_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_sortedvelocity_log_runSort_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/03_velocity_log_runSort_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.SORT, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_SORTPOSTPASS() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_density_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_gridcellindex_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_gridcellindexfixedup_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_index_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_indexback_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_neighbormap_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_position_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_pressure_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_sortedposition_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_sortedvelocity_log_runSortPostPass_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/04_velocity_log_runSortPostPass_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.SORT_POST_PASS, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_RUNINDEX() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_density_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_gridcellindex_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_gridcellindexfixedup_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_index_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_indexback_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_neighbormap_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_position_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_pressure_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_sortedposition_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_sortedvelocity_log_runIndexx_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/05_velocity_log_runIndexx_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.INDEX, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_RUNINDEXPOSTPASS() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_density_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_gridcellindex_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_gridcellindexfixedup_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_index_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_indexback_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_neighbormap_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_position_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_pressure_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_sortedposition_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_sortedvelocity_log_runIndexPostPass_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/06_velocity_log_runIndexPostPass_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.INDEX_POST_PASS, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_780_FINDNEIGHBOR() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_density_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_gridcellindex_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_gridcellindexfixedup_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_index_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_indexback_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_neighbormap_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_position_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_pressure_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_sortedposition_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_sortedvelocity_log_runFindNeighbors_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/liquid_780/checkpoints/step1/07_velocity_log_runFindNeighbors_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.FIND_NEIGHBORS, logs, this.getClass().getResource("/sphModel_liquid_780.xml"));
	}
	
	@Test
	public void testCheckpoints_1575_CLEARBUFFERS() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_density_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_gridcellindex_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_gridcellindexfixedup_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_index_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_indexback_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_neighbormap_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_position_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_pressure_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_sortedposition_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_sortedvelocity_log_runClearBuffers_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/01_velocity_log_runClearBuffers_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.CLEAR_BUFFERS, logs, this.getClass().getResource("/sphModel_elastic_1575.xml"));
	}
	
	@Test
	public void testCheckpoints_1575_HASHPARTICLES() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();
		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_density_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_gridcellindex_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_gridcellindexfixedup_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_index_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_indexback_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_neighbormap_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_position_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_pressure_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_sortedposition_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_sortedvelocity_log_runHashParticles_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/02_velocity_log_runHashParticles_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.HASH_PARTICLES, logs, this.getClass().getResource("/sphModel_elastic_1575.xml"));
	}
	
	@Test
	public void testCheckpoints_1575_COMPUTEELASTIC() throws Exception {
		// load reference values at various steps from C++ version
		Map<BuffersEnum, URL> logs = new LinkedHashMap<BuffersEnum, URL>();

		logs.put(BuffersEnum.RHO, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_density_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_gridcellindex_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.GRID_CELL_INDEX_FIXED, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_gridcellindexfixedup_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_index_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.PARTICLE_INDEX_BACK, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_indexback_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.NEIGHBOR_MAP, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_neighbormap_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_position_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.PRESSURE, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_pressure_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.SORTED_POSITION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_sortedposition_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.SORTED_VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_sortedvelocity_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.VELOCITY, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_velocity_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.ACCELERATION, StepValidationTest.class.getResource("/results/elastic_1575/checkpoints/step1/10_acceleration_log_run_pcisph_computeElasticForces_0.txt"));
		logs.put(BuffersEnum.ELASTIC_CONNECTIONS, StepValidationTest.class.getResource("/results/elastic_1575/elastic_log_0.txt"));
		
		evaluateCheckpoint(KernelsEnum.COMPUTE_ELASTIC_FORCES, logs, this.getClass().getResource("/sphModel_elastic_1575.xml"));
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
	
	/*
	 * Retrieves an array of ints from a tab separated string of values
	 * */
	private Integer[] getIntValues(String values)
	{
		String[] series = values.split("\t");
		Integer[] intSeries = new Integer[series.length];
		
		for(int i=0; i<series.length; i++)
		{
			// equivalent to -1 in C++
			if(series[i].equals("4294967295")) series[i] = "-1";
			intSeries[i] = Integer.parseInt(series[i]);
		}
		
		return intSeries;
	}
	
	/*
	 * Retrieves an array of floats from a tab separated string of values
	 * */
	private Float[] getFloatValues(String values)
	{
		String[] series = values.split("\t");
		Float[] floatSeries = new Float[series.length];
		
		for(int i=0; i<series.length; i++)
		{
			floatSeries[i] = Float.parseFloat(series[i]);
		}
		
		return floatSeries;
	}
	
	/*
	 * Evaluates checkpoints
	 * NOTE: this routine knows a whole about log file formats and stuff - it's ugly but it works
	 * */
	private void evaluateCheckpoint(KernelsEnum checkpoint, Map<BuffersEnum, URL> logs, URL modelURL) throws Exception
	{
		// load reference values at various steps from C++ version
		Map<BuffersEnum, String[]> checkpointReferenceValuesMap = new HashMap<BuffersEnum, String[]>();
		for (Map.Entry<BuffersEnum, URL> entry : logs.entrySet())
		{
			String fileContent = PCISPHTestUtilities.readFile(logs.get(entry.getKey()).getPath());
			checkpointReferenceValuesMap.put(entry.getKey(), fileContent.split(System.getProperty("line.separator")));
		}
		
		// load Java generated scene
		SPHModelInterpreterService modelInterpreter = new SPHModelInterpreterService();
		SPHModelX model = (SPHModelX)modelInterpreter.readModel(modelURL,null,"");
		
		SPHSolverService solver = new SPHSolverService(true);
		solver.initialize(model);
		
		Map<BuffersEnum, Integer> mismatchingValuesPerBuffers = new LinkedHashMap<BuffersEnum, Integer>();

		AspectNode aspect= new AspectNode("Aspect");
		// calculate step
		solver.solve(new TimeConfiguration(null, 1, null),aspect);
		
		// get checkpoint of interest
		PCISPHCheckPoint checkpoint_values = solver.getCheckpointsMap().get(checkpoint);
		
		// get buffer sizes
		Map<BuffersEnum, Integer> dimensions = solver.getBuffersSizeMap();
		
		// compare buffer values at first checkpoint
		for (Map.Entry<BuffersEnum, String[]> entry : checkpointReferenceValuesMap.entrySet())
		{
			switch (entry.getKey()) {
		        case RHO:  
		        	List<Float> rho_calculatedValues = checkpoint_values.rho;
		        	int rho_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length *2);
		        	int j = 0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 2)
		        	{	
		        		// get vector from ref values
		        		Vector3D vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(vector.getX().floatValue() != rho_calculatedValues.get(i).floatValue()) 
		        		{
		        			rho_mismatches++;
		        		}
		        		if(vector.getY().floatValue() != rho_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			rho_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), rho_mismatches);
		            break;
		        case ELASTIC_CONNECTIONS:  
		        	// liquid scene has no elastic connections - log files are empty
		        	List<Float> el_calculatedValues = checkpoint_values.elasticConnections;
		        	int el_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 4);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 4)
		        	{	
		        		// get vector from ref values
		        		Vector3D el_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(el_vector.getX().floatValue() != el_calculatedValues.get(i).floatValue()) 
		        		{
		        			el_mismatches++;
		        		}
		        		if(el_vector.getY().floatValue() != el_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			el_mismatches++;
		        		}
		        		if(el_vector.getZ().floatValue() != el_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			el_mismatches++;
		        		}
		        		if(el_vector.getP().floatValue() != el_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			el_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), el_mismatches);
		            break;
		        case GRID_CELL_INDEX:  
		        	List<Integer> grid_cell_calculatedValues = checkpoint_values.gridCellIndex;
	        		
		        	// get array of values from ref values - this buffer is all on one line
	        		Integer[] grid_cell_vector = getIntValues(checkpointReferenceValuesMap.get(entry.getKey())[0]);
		        	
	        		int grid_cell_fixed_mismatches = 0;
	        		Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey())[0].split("\t").length);
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i++)
		        	{	
		        		if(grid_cell_vector[i].intValue() != grid_cell_calculatedValues.get(i).intValue()) 
		        		{
		        			grid_cell_fixed_mismatches++;
		        		}
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), grid_cell_fixed_mismatches);
		            break;
		        case GRID_CELL_INDEX_FIXED:  
		        	List<Integer> grid_idx_calculatedValues = checkpoint_values.gridCellIndexFixedUp;
	        		
		        	// get array of values from ref values - this buffer is all on one line
	        		Integer[] vector = getIntValues(checkpointReferenceValuesMap.get(entry.getKey())[0]);
		        	
	        		int grid_idx_fixed_mismatches = 0;
	        		Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey())[0].split("\t").length);
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i++)
		        	{	
		        		if(vector[i].intValue() != grid_idx_calculatedValues.get(i).intValue()) 
		        		{
		        			grid_idx_fixed_mismatches++;
		        		}
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), grid_idx_fixed_mismatches);
		            break;
		        case PARTICLE_INDEX:  
		        	List<Integer> index_calculatedValues = checkpoint_values.particleIndex;
	        		
		        	// get array of values from ref values - this buffer is all on one line
	        		Integer[] index_vector = getIntValues(checkpointReferenceValuesMap.get(entry.getKey())[0]);
		        	
	        		int index_mismatches = 0;
	        		Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey())[0].split("\t").length);
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i++)
		        	{	
		        		if(index_vector[i].intValue() != index_calculatedValues.get(i).intValue()) 
		        		{
		        			index_mismatches++;
		        		}
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), index_mismatches);
		            break;
		        case PARTICLE_INDEX_BACK:  
		        	List<Integer> index_back_calculatedValues = checkpoint_values.particleIndexBack;
	        		
		        	// get array of values from ref values - this buffer is all on one line
	        		Integer[] index_back_vector = getIntValues(checkpointReferenceValuesMap.get(entry.getKey())[0]);
		        	
	        		int index_back_mismatches = 0;
	        		Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey())[0].split("\t").length);
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i++)
		        	{	
		        		if(index_back_vector[i].intValue() != index_back_calculatedValues.get(i).intValue()) 
		        		{
		        			index_back_mismatches++;
		        		}
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), index_back_mismatches);
		            break;
		        case NEIGHBOR_MAP:  
		        	List<Float> nm_calculatedValues = checkpoint_values.neighborMap;
		        	int nm_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length *2);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 2)
		        	{	
		        		// get vector from ref values
		        		Vector3D nm_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(nm_vector.getX().floatValue() != nm_calculatedValues.get(i).floatValue()) 
		        		{
		        			nm_mismatches++;
		        		}
		        		if(nm_vector.getY().floatValue() != nm_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			nm_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), nm_mismatches);
		            break;
		        case POSITION:  
		        	List<Float> pos_calculatedValues = checkpoint_values.position;
		        	int pos_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 4);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 4)
		        	{	
		        		// get vector from ref values
		        		Vector3D pos_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(pos_vector.getX().floatValue() != pos_calculatedValues.get(i).floatValue()) 
		        		{
		        			pos_mismatches++;
		        		}
		        		if(pos_vector.getY().floatValue() != pos_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			pos_mismatches++;
		        		}
		        		if(pos_vector.getZ().floatValue() != pos_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			pos_mismatches++;
		        		}
		        		if(pos_vector.getP().floatValue() != pos_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			pos_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), pos_mismatches);
		            break;
		        case ACCELERATION:  
		        	List<Float> acc_calculatedValues = checkpoint_values.acceleration;
		        	int acc_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 8);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 8)
		        	{	
		        		// get vector from ref values
		        		Float[] acc_vector = getFloatValues(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(acc_vector[0].floatValue() != acc_calculatedValues.get(i).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[1].floatValue() != acc_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[2].floatValue() != acc_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[3].floatValue() != acc_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[4].floatValue() != acc_calculatedValues.get(i + 4).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[5].floatValue() != acc_calculatedValues.get(i + 5).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[6].floatValue() != acc_calculatedValues.get(i + 6).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		if(acc_vector[7].floatValue() != acc_calculatedValues.get(i + 7).floatValue()) 
		        		{
		        			acc_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), acc_mismatches);
		            break;
		        case PRESSURE:  
		        	List<Float> press_calculatedValues = checkpoint_values.pressure;
		        	int press_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 4);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 4)
		        	{	
		        		// get vector from ref values
		        		Vector3D press_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(press_vector.getX().floatValue() != press_calculatedValues.get(i).floatValue()) 
		        		{
		        			press_mismatches++;
		        		}
		        		if(press_vector.getY().floatValue() != press_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			press_mismatches++;
		        		}
		        		if(press_vector.getZ().floatValue() != press_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			press_mismatches++;
		        		}
		        		if(press_vector.getP().floatValue() != press_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			press_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), press_mismatches);
		            break;
		        case SORTED_POSITION:  
		        	List<Float> sort_pos_calculatedValues = checkpoint_values.sortedPosition;
		        	int sort_pos_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 8);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 8)
		        	{	
		        		// get vector from ref values
		        		Float[] sort_pos_vector = getFloatValues(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(sort_pos_vector[0].floatValue() != sort_pos_calculatedValues.get(i).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[1].floatValue() != sort_pos_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[2].floatValue() != sort_pos_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[3].floatValue() != sort_pos_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[4].floatValue() != sort_pos_calculatedValues.get(i + 4).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[5].floatValue() != sort_pos_calculatedValues.get(i + 5).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[6].floatValue() != sort_pos_calculatedValues.get(i + 6).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		if(sort_pos_vector[7].floatValue() != sort_pos_calculatedValues.get(i + 7).floatValue()) 
		        		{
		        			sort_pos_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), sort_pos_mismatches);
		            break;
		        case SORTED_VELOCITY:  
		        	List<Float> sort_vel_calculatedValues = checkpoint_values.sortedVelocity;
		        	int sort_vel_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 4);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 4)
		        	{	
		        		// get vector from ref values
		        		Vector3D sort_vel_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(sort_vel_vector.getX().floatValue() != sort_vel_calculatedValues.get(i).floatValue()) 
		        		{
		        			sort_vel_mismatches++;
		        		}
		        		if(sort_vel_vector.getY().floatValue() != sort_vel_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			sort_vel_mismatches++;
		        		}
		        		if(sort_vel_vector.getZ().floatValue() != sort_vel_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			sort_vel_mismatches++;
		        		}
		        		if(sort_vel_vector.getP().floatValue() != sort_vel_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			sort_vel_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), sort_vel_mismatches);
		            break;
		        case VELOCITY:  
		        	List<Float> vel_calculatedValues = checkpoint_values.velocity;
		        	int vel_mismatches = 0;
		        	Assert.assertTrue(dimensions.get(entry.getKey()).intValue() == checkpointReferenceValuesMap.get(entry.getKey()).length * 4);
		        	j=0;
		        	for(int i = 0; i < dimensions.get(entry.getKey()); i = i + 4)
		        	{	
		        		// get vector from ref values
		        		Vector3D vel_vector = get3DVector(checkpointReferenceValuesMap.get(entry.getKey())[j]);
		        		
		        		// it sucks, but all the if statements are separate to facilitate debugging
		        		if(vel_vector.getX().floatValue() != vel_calculatedValues.get(i).floatValue()) 
		        		{
		        			vel_mismatches++;
		        		}
		        		if(vel_vector.getY().floatValue() != vel_calculatedValues.get(i + 1).floatValue()) 
		        		{
		        			vel_mismatches++;
		        		}
		        		if(vel_vector.getZ().floatValue() != vel_calculatedValues.get(i + 2).floatValue()) 
		        		{
		        			vel_mismatches++;
		        		}
		        		if(vel_vector.getP().floatValue() != vel_calculatedValues.get(i + 3).floatValue()) 
		        		{
		        			vel_mismatches++;
		        		}
		        		
		        		j++;
		        	}
		        	// record mismatch
        			mismatchingValuesPerBuffers.put(entry.getKey(), vel_mismatches);
		            break;
		        default:
		        	// do nothing
		        	break;
			}
		}
		
		
		StringBuilder sb = new StringBuilder();
		// compare buffer values at first checkpoint
		for (Map.Entry<BuffersEnum, Integer> entry : mismatchingValuesPerBuffers.entrySet())
		{
			BuffersEnum buffer = entry.getKey();
			// assert there are no mismatching values
			if (entry.getValue() != 0)
				sb.append(entry.getValue()).append(" of ").append(dimensions.get(buffer))
				   .append(" mismatching values found on buffer ").append(buffer.toString()).append("\n");
		}
		
		String msg = sb.toString();
		if (!msg.isEmpty())
			Assert.fail(msg);
	}
}
