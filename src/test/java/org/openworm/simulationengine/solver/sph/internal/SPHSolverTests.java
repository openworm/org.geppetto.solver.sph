package org.openworm.simulationengine.solver.sph.internal;

import static java.lang.System.out;

import java.io.IOException;
import org.bridj.Pointer;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import org.openworm.simulationengine.core.constants.PhysicsConstants;
import org.openworm.simulationengine.solver.sph.SPHContstants;
import org.openworm.simulationengine.solver.sph.SPHSolverService;



public class SPHSolverTests {

	private static SPHSolverService solver = new SPHSolverService();


	static final int NO_PARTICLE_ID = -1;
	static final float EPSILON = 0.1f;
	static final int MAX_ERRORS = 10;

	@BeforeClass 
	public static void InitSover() {
		try {
			solver.init();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	@AfterClass 
	public static void ReleaseSolver(){
		solver.cleanContext();
	}
	

	@Test
	public void _runStep(){
		testClearBuffers();
		testHashParticles();
		testSort();
		testSortPostPass();		
		testIndexx();
		testIndexPostPass();
		testFindNeighbors();
		testComputeDensityPressure();
		testComputeAcceleration();
		testIntegrate();
	}
	
	//@Test
	public void testClearBuffers(){
		solver.queue.finish();
		System.out.println("Testing stage ClearBuffers:");
		boolean thisTestPassed = true;
		solver._runClearBuffers();
		solver.queue.finish();
		thisTestPassed = testEndClearBuffers();
		if (!thisTestPassed)
			System.out.println("Test FAILED in stage ClearBuffers");
		else
			System.out.println("Test SUCCEEDED in stage ClearBuffers");
	}
	
	//@Test
	public void testHashParticles(){
		// Invoke stage HashParticles
		solver.queue.finish();
		boolean thisTestPassed = true;
		solver._runHashParticles();
		solver.queue.finish();
		thisTestPassed = testEndHashParticles();

		if( !thisTestPassed ){
			System.out.println("Test FAILED in stage HashParticles");
		}else{
			System.out.println("Test SUCCEEDED in stage HashParticles");
		}
		
	}
	
	//@Test
	public void testSort(){
		// Invoke stage Sort
		boolean thisTestPassed = true;
		out.println("Testing stage Sort:");
		thisTestPassed = testBeginSort();
		solver._runSort();
		thisTestPassed = testEndSort();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage Sort");
		}else{
			out.println("Test SUCCEEDED in stage Sort" );
		}
	}
	
	//@Test
	public void testSortPostPass(){
		// Invoke stage SortPostPass
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage SortPostPass:");
		solver._runSortPostPass();
		solver.queue.finish();
		thisTestPassed = testEndSortPostPass();

		if( !thisTestPassed ){
			out.println("Test FAILED in stage SortPostPass");
		}else{
			out.println("Test SUCCEEDED in stage SortPostPass");
		}
	}
	
	//@Test
	public void testIndexx(){
		// Invoke stage Indexx
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage Indexx:");
		solver._runIndexx();
		solver.queue.finish();
		thisTestPassed = testEndIndexx();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage Indexx");
		}else{
			out.println("Test SUCCEEDED in stage Indexx");
		}
	}
	
	//@Test
	public void testIndexPostPass(){
		// Invoke stage IndexPostPass
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage IndexPostPass:");
		solver._runIndexPostPass();
		solver.queue.finish();
		thisTestPassed = testEndIndexPostPass();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage IndexPostPass");
		}else{
			out.println("Test SUCCEEDED in stage IndexPostPass");
		}
	}
	
	//@Test
	public void testFindNeighbors(){
		// Invoke stage FindNeighbors
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage FindNeighbors:");
		solver._runFindNeighbors();
		solver.queue.finish();
		thisTestPassed = testEndFindNeighbors();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage FindNeighbors");
		}else{
			out.println("Test SUCCEEDED in stage FindNeighbors");
		}
	}
	
	//@Test
	public void testComputeDensityPressure(){
		// Invoke stage ComputeDensityPressure
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage ComputeDensityPressure:");
		solver._runComputeDensityPressure();
		solver.queue.finish();
		thisTestPassed = testEndComputeDensityPressure();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage ComputeDensityPressure");
		}else{
			out.println("Test SUCCEEDED in stage ComputeDensityPressure");
		}
	}
	
	//@Test
	public void testComputeAcceleration(){
		// Invoke stage ComputeAcceleration
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage ComputeAcceleration:");
		solver._runComputeAcceleration();
		solver.queue.finish();
		thisTestPassed = testEndComputeAcceleration();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage ComputeAcceleration");
		}else{
			out.println("Test SUCCEEDED in stage ComputeAcceleration");
		}
	}
	
	//@Test
	public void testIntegrate(){
		// Invoke stage Integrate
		solver.queue.finish();
		boolean thisTestPassed = true;
		out.println("Testing stage Integrate:");
		solver._runIntegrate();
		solver.queue.finish();
		thisTestPassed = testEndIntegrate();
		if( !thisTestPassed ){
			out.println("Test FAILED in stage Integrate");
		}else{
			out.println("Test SUCCEEDED in stage Integrate");
		}
	}
	
	private boolean testEndClearBuffers(){
		//test if all ellement of neighborMap is -1
		boolean result = true;
		int messageCount = 0;
		solver.neighborMapPtr = solver.neighborMap.read(solver.queue);
		solver.queue.finish();
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			float neighborId = solver.neighborMapPtr.get(i*2);//[ i * 2 ];
			if( !( neighborId == -1 )){
				System.out.println("Error: failed to initialize neighborMap at location " + i);
				result = false;
			}
			if(result == false )
				messageCount += 1;
			if( messageCount > MAX_ERRORS ) break;
		}//for
		
		return result;
	}
	
	boolean testEndHashParticles(){
		boolean result = true;
		int messageCount = 0;
		solver.position.read(solver.queue, solver.positionPtr, true);
		solver.queue.finish();
		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = (int)solver.particleIndexPtr.get(i * 2);
			int particleId = (int)solver.particleIndexPtr.get(i * 2 + 1);
			if ( cellId < 0 || cellId >= solver.gridCellCount ){
				System.out.println("Error: bad cellId " + cellId + " at particle " + particleId);
				result = false;
				messageCount++;
			}
			int ix = (int)(( solver.positionPtr.get(4 * particleId) - SPHContstants.XMIN_F ) / PhysicsConstants.HASH_GRID_CELL_SIZE ) % solver.gridCellsX;
			int iy = (int)(( solver.positionPtr.get(4 * particleId + 1) - SPHContstants.YMIN_F ) / PhysicsConstants.HASH_GRID_CELL_SIZE ) % solver.gridCellsY;
			int iz = (int)(( solver.positionPtr.get(4 * particleId + 2) - SPHContstants.ZMIN_F ) / PhysicsConstants.HASH_GRID_CELL_SIZE ) % solver.gridCellsZ;
			if( ix < 0 ) ix += solver.gridCellsX;
			if( iy < 0 ) iy += solver.gridCellsY;
			if( iz < 0 ) iz += solver.gridCellsZ;
			assert( ix >= 0 && iy >= 0 && iz >= 0 );
			int computedCellId = ix + iy * solver.gridCellsX + iz * solver.gridCellsX * solver.gridCellsY;
			computedCellId = computedCellId & 0xffff; // truncate to low 16 bits
			if( computedCellId != cellId ){
				System.out.println("Error: inconsistent cellId "+ cellId +" at particle " +particleId+" expected " + computedCellId);
				result = false;
				messageCount++;
			}
		}
		return result;
	}
	
	
	boolean testBeginSort(){
		int messageCount = 0;
		boolean result = true;

		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = solver.particleIndexPtr.get(i * 2);
			int particleId = solver.particleIndexPtr.get(i * 2 + 1);
			if( cellId < 0 || cellId > solver.gridCellCount ){
				System.err.println("Error: bad cellId " + cellId + " input data to sort, element " + i);
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
			}
			if( particleId != i ){
				System.err.println("Error: bad particleId "+ particleId +" input data to sort, element "+i);
				result = false;
				if( ++messageCount > MAX_ERRORS ) break;
				continue;
			}
		}//for
		return result;
	}
	
	private boolean testEndSort() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		int[] test = new int[SPHContstants.PARTICLE_COUNT * 2];
		for(int i=0;i < SPHContstants.PARTICLE_COUNT *2;i++){
			test[i] = solver.particleIndexPtr.get(i);
		}
		int lastCellId = 0;
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = solver.particleIndexPtr.get(i * 2);
			int particleId = solver.particleIndexPtr.get(i * 2 + 1);
			if ( cellId < lastCellId || cellId < 0 ){
				out.println("Error: data out of order "+cellId+" at position "+i);
				result = false;
			}
			lastCellId = cellId;
			if(result == false)
				messageCount += 1;
			if( messageCount > MAX_ERRORS ) break;
		}//for
		return result;
	}
	
	private boolean testEndSortPostPass() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		solver.position.read( solver.queue, solver.positionPtr, true);
		solver.queue.finish();
		solver.velocity.read(solver.queue, solver.velocityPtr, true);
		solver.queue.finish();
		solver.sortedPosition.read(solver.queue, solver.sortedPositionPtr, true);
		solver.queue.finish();
		solver.sortedVelocity.read(solver.queue, solver.sortedVelocityPtr, true);
		solver.queue.finish();

		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = solver.particleIndexPtr.get(2*i);
			int serialId = solver.particleIndexPtr.get(2*i + 1);
			if( cellId < 0 || cellId > solver.gridCellCount
				|| serialId < 0 || serialId > SPHContstants.PARTICLE_COUNT ){
					out.println("Error: impossible results from sort, particle "+i 
						+" cell "+cellId+" serialId "+serialId);
					result = false;
					if( ++messageCount > MAX_ERRORS ) break;
					continue;
			}
			float[] oldPosition = readField( 4 * serialId, solver.positionPtr );
			float[] newPosition = readField(4 * i, solver.sortedPositionPtr );
			if( oldPosition[ 0 ] != newPosition[ 0 ] ||
				oldPosition[ 1 ] != newPosition[ 1 ] ||
				oldPosition[ 2 ] != newPosition[ 2 ] ){
					out.println("Error: particle "+ i +" position does not match sorting element "
						+ serialId);
					result = false;
					if( ++messageCount > MAX_ERRORS ) break;
					continue;
			}
			float[] oldVelocity = readField( 4 * serialId, solver.velocityPtr );
			float[] newVelocity = readField( 4 * i, solver.sortedVelocityPtr );
			if( oldVelocity[ 0 ] != newVelocity[ 0 ] ||
				oldVelocity[ 1 ] != newVelocity[ 1 ] ||
				oldVelocity[ 2 ] != newVelocity[ 2 ] ){
					out.println("Error: particle "+ i + " velocity does not match sorting element "
						+ serialId);
					result = false;
					if( ++messageCount > MAX_ERRORS ) break;
					continue;
			}
		}//for
		return result;
	}
	
	boolean testEndComputeAcceleration( ){

		boolean result = true;
		int messageCount = 0;
		
		solver.neighborMap.read(solver.queue, solver.neighborMapPtr, true);
		solver.queue.finish();
		solver.pressure.read(solver.queue, solver.pressurePtr, true);
		solver.queue.finish();
		solver.rho.read(solver.queue, solver.rhoPtr, true);
		solver.queue.finish();
		solver.rhoInv.read(solver.queue, solver.rhoInvPtr, true);
		solver.queue.finish();
		solver.sortedPosition.read(solver.queue, solver.sortedPositionPtr, true);
		solver.queue.finish();
		solver.sortedVelocity.read(solver.queue, solver.sortedVelocityPtr, true);
		solver.queue.finish();
		solver.acceleration.read(solver.queue, solver.accelerationPtr, true);
		solver.queue.finish();
		
		float[] gravity = { PhysicsConstants.GRAVITY_X, PhysicsConstants.GRAVITY_Y, PhysicsConstants.GRAVITY_Z };
		float totalMagnitude = 0.0f;
		float hScaled = PhysicsConstants.H * PhysicsConstants.SIMULATION_SCALE;
		
		for( int particleId = 0; particleId < SPHContstants.PARTICLE_COUNT; ++particleId ){
			float rho_i = solver.rhoPtr.get(particleId);
			float rho_i_inv = 1.0f / rho_i;
			float p_i = solver.pressurePtr.get(particleId);
			//float * nm = _neighborMap + particleId * NEIGHBOR_COUNT * 2;
			float[] computedAcceleration = { 0.0f, 0.0f, 0.0f }; 
			float[] position_i = readField( particleId * 4, solver.sortedPositionPtr );
			float[] v_i = readField( particleId * 4, solver.sortedVelocityPtr );
			float[] gradP = { 0.0f, 0.0f, 0.0f };
			float[] del2V = { 0.0f, 0.0f, 0.0f };
		
		// integrate the gradP and del^2 V terms over all particles
		
			for( int neighborNum = 0 /* -1 */; neighborNum < SPHContstants.NEIGHBOR_COUNT; ++neighborNum ){
				int neighborParticleId;
				float r;
				float rho_j_inv;
				float p_j;
				float[] neighborPosition;
				float[] d = new float[3];
				float[] v_j;
				float t = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + neighborNum * 2);
				neighborParticleId = (int)t;
				if( neighborParticleId == NO_PARTICLE_ID ) continue;
				r = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + neighborNum * 2 + 1);
				float rhoP_j = solver.rhoPtr.get(neighborParticleId);
				rho_j_inv = 1.0f / rhoP_j;
				p_j = solver.pressurePtr.get(neighborParticleId);
				neighborPosition = readField( neighborParticleId * 4, solver.sortedPositionPtr );
				d[ 0 ] = position_i[ 0 ] - neighborPosition[ 0 ];
				d[ 1 ] = position_i[ 1 ] - neighborPosition[ 1 ];
				d[ 2 ] = position_i[ 2 ] - neighborPosition[ 2 ];
				d[ 0 ] *= PhysicsConstants.SIMULATION_SCALE;
				d[ 1 ] *= PhysicsConstants.SIMULATION_SCALE;
				d[ 2 ] *= PhysicsConstants.SIMULATION_SCALE;
				d[ 0 ] /= r;
				d[ 1 ] /= r;
				d[ 2 ] /= r;
				v_j =  readField( neighborParticleId * 4, solver.sortedVelocityPtr );		  
		
				// gradP
				
				float x = hScaled - r;
				float[] gradWspiky = new float[ 3 ];
				gradWspiky[ 0 ] = x * x * d[ 0 ] * PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT;
				gradWspiky[ 1 ] = x * x * d[ 1 ] * PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT;
				gradWspiky[ 2 ] = x * x * d[ 2 ] * PhysicsConstants.GRAD_W_SPIKY_COEFFICIENT;
				for( int j = 0; j < 3; ++j ){
					gradP[ j ] += PhysicsConstants.MASS * ( p_i + p_j ) * 0.5f * rho_j_inv * gradWspiky[ j ];
				}//for
				
				// del^2 V
				
				float del2Wviscosity = ( hScaled - r ) * PhysicsConstants.DEL_2_W_VISCOSITY_COEFFICIENT;
				d[ 0 ] = ( v_j[ 0 ] - v_i[ 0 ] );
				d[ 1 ] = ( v_j[ 1 ] - v_i[ 1 ] );
				d[ 2 ] = ( v_j[ 2 ] - v_i[ 2 ] );
				for( int j = 0; j < 3; ++j ){
					del2V[ j ] += PhysicsConstants.MASS * d[ j ] * rho_j_inv * del2Wviscosity;
				}//for
			}//for
		
			for( int j = 0; j < 3; ++j ){
				computedAcceleration[ j ] = rho_i_inv * ( PhysicsConstants.MU * del2V[ j ] - gradP[ j ] );
			}//for
			// apply CFL limiting
			float magnitudeSquared = computedAcceleration[ 0 ] * computedAcceleration[ 0 ]
					+ computedAcceleration[ 1 ] * computedAcceleration[ 1 ]
							+ computedAcceleration[ 2 ] * computedAcceleration[ 2 ];
			float magnitude = (float)Math.sqrt( magnitudeSquared );
			if( magnitude > PhysicsConstants.CFLLimit ){
				float scale = PhysicsConstants.CFLLimit / magnitude;
				for( int j = 0; j < 3; ++j ){
					computedAcceleration[ j ] *= scale;
				}//for
			}
			totalMagnitude += magnitude;
		
			float[] acceleration = readField( particleId * 4, solver.accelerationPtr );
			float[] d = new float[ 3 ];
			d[ 0 ] = computedAcceleration[ 0 ] - acceleration[ 0 ];
			d[ 1 ] = computedAcceleration[ 1 ] - acceleration[ 1 ];
			d[ 2 ] = computedAcceleration[ 2 ] - acceleration[ 2 ];
			float l2Norm = (float)Math.sqrt( d[ 0 ] * d[ 0 ] + d[ 1 ] * d[ 1 ] + d[ 2 ] * d[ 2 ] );
			if( l2Norm > EPSILON ){
				System.err.println("Error: particle "+particleId+" has unexpected acceleration ("
						+acceleration[ 0 ]+ ", "+ acceleration[ 1 ]+ ", "+acceleration[ 2 ]+") \texpected ("+ computedAcceleration[ 0 ]+", "
						+computedAcceleration[ 1 ]+ ", "+computedAcceleration[ 2 ]+")");
				result = false;
				messageCount++;
			}
			if( messageCount > MAX_ERRORS ) break;
		}//for
		
		float[] totalAcceleration = { 0.0f, 0.0f, 0.0f };
		
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			float[] a = readField(i*4, solver.accelerationPtr);
			totalAcceleration[ 0 ] += a[ 0 ];
			totalAcceleration[ 1 ] += a[ 1 ];
			totalAcceleration[ 2 ] += a[ 2 ];
		}//for
		
		float[] meanAcceleration = new float[ 3 ];
		for( int j = 0; j < 3; ++j ){
			meanAcceleration[ j ] = totalAcceleration[ j ] / SPHContstants.PARTICLE_COUNT;
		}//for
		float magnitude = (float)Math.sqrt( 
				meanAcceleration[ 0 ] * meanAcceleration[ 0 ] +
				meanAcceleration[ 1 ] * meanAcceleration[ 1 ] +
				meanAcceleration[ 2 ] * meanAcceleration[ 2 ]
				);
		float meanMagnitude = totalMagnitude / SPHContstants.PARTICLE_COUNT;
		System.err.println("mean acceleration ("+meanAcceleration[ 0 ]+","+ meanAcceleration[ 1 ]+","+meanAcceleration[ 2 ]+") mean magnitude "+ meanMagnitude);
		return result;
	}
	
	
	private boolean testEndComputeDensityPressure() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.neighborMap.read(solver.queue,solver.neighborMapPtr, true);
		solver.queue.finish();
		solver.pressure.read(solver.queue, solver.pressurePtr, true);
		solver.queue.finish();
		solver.rho.read(solver.queue, solver.rhoPtr, true);
		solver.queue.finish();
		solver.rhoInv.read(solver.queue, solver.rhoInvPtr, true);
		solver.queue.finish();


		float hScaled = PhysicsConstants.H * PhysicsConstants.SIMULATION_SCALE;

		for( int particleId = 0; particleId < SPHContstants.PARTICLE_COUNT; ++particleId ){
			//float[] nm = readField(particleId * SPHSolverService.NEIGHBOR_COUNT * 2, solver.neighborMap) ;
			float computedDensity = 0.0f;

			for( int neighborNum = 0; neighborNum < SPHContstants.NEIGHBOR_COUNT; ++neighborNum ){
				float r = 0.0f;
				float neighborParticleId = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + neighborNum * 2);//nm[0];
				if( neighborParticleId == NO_PARTICLE_ID ) continue;
				float distance = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + neighborNum * 2 + 1);
				r = distance;
				float x = hScaled * hScaled - r * r;
				float Wpoly6 = x * x * x * PhysicsConstants.W_POLY_6_COEFFICIENT;
				float contribution = PhysicsConstants.MASS * Wpoly6;
				computedDensity += contribution;
			}//for

			float drho = computedDensity - PhysicsConstants.RHO0;
			float k = PhysicsConstants.STIFFNESS;
			float computedPressure = k * drho;

			float density = solver.rhoPtr.get(particleId);
			float p = solver.pressurePtr.get(particleId);

			float densityError = 1.0f - density / computedDensity;
			//float densityError = fabs( density - computedDensity );
			if( densityError > EPSILON ){
				out.println("Error: particle "+particleId+" expected density "
					+computedDensity+" but kernel computed "+density+" error "+densityError+"%");
				result = false;
				messageCount++;
			}
			float pressureError = 1.0f - p / computedPressure;
			//float pressureError = fabs( p - computedPressure );
			if( pressureError > EPSILON ){
				out.println("Error: particle "+particleId+" expected pressure "
					+computedPressure+" but kernel computed "+p+" error "+pressureError+"%");
				result = false;
				messageCount++;
			}
			if( messageCount > MAX_ERRORS ) break;
		}//for

		float totalRho = 0.0f;
		float totalP = 0.0f;

		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			totalRho += solver.rhoPtr.get(i);
			totalP += solver.pressurePtr.get(i);
		}//for

		float meanRho = totalRho / SPHContstants.PARTICLE_COUNT;
		float meanP = totalP / SPHContstants.PARTICLE_COUNT;
		out.println("mean rho "+meanRho+" mean pressure "+meanP);
		return result;
	}
		
	private boolean testEndFindNeighbors() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.gridCellIndexFixedUp.read( solver.queue, solver.gridCellIndexFixedUpPtr, true);
		solver.sortedPosition.read(solver.queue, solver.sortedPositionPtr, true);
		solver.neighborMap.read( solver.queue, solver.neighborMapPtr, true);
		solver.queue.finish();

		// make sure neighbors are symmetric and agree with distance

		for( int particleId = 0; particleId < SPHContstants.PARTICLE_COUNT; ++particleId ){
			float[] nm = readField(particleId * SPHContstants.NEIGHBOR_COUNT * 2, solver.neighborMapPtr);
			float[] myPosition = readField( 4*particleId, solver.sortedPositionPtr );

			for( int neighborNum = 0; neighborNum < SPHContstants.NEIGHBOR_COUNT; ++neighborNum ){
				float neighborParticleId = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + neighborNum * 2);
				if( (int)neighborParticleId == NO_PARTICLE_ID ) continue;

				// make sure this particle is listed in the neighbors neighbor map
				boolean isNeighbor = false;
				for( int j = 0; j < SPHContstants.NEIGHBOR_COUNT; ++j ){
					int nnNum = (int)((float)solver.neighborMapPtr.get(j * 2 ));
					isNeighbor |= nnNum == particleId;
				}//for
				if( !isNeighbor ){
					out.println("Warning: particle "+particleId+" has neighbor particle "
						+neighborParticleId+" but is not listed in the neighbor map for that particle\n");
					messageCount++;
					continue;
				}

				// make sure the computed distance is correct

				float distance = solver.neighborMapPtr.get(neighborNum * 2 + 1);
				float[] neighborPosition = readField( ((int)neighborParticleId) *4, solver.sortedPositionPtr );

				float dx = myPosition[ 0 ] - neighborPosition[ 0 ];
				float dy = myPosition[ 1 ] - neighborPosition[ 1 ];
				float dz = myPosition[ 2 ] - neighborPosition[ 2 ];
				float distanceSquared = dx * dx + dy * dy + dz * dz;
				float computedDistance = (float)Math.sqrt( distanceSquared ) * PhysicsConstants.SIMULATION_SCALE;

				if( Math.abs( 1.0f - computedDistance / distance ) > EPSILON ){
					out.println("Error: particle "+particleId+" has neighbor particle "
						+neighborParticleId+" but distance "+distance
						+" does not agree with computed distance "
						+computedDistance);
					result = false;
					messageCount++;
				}	
				if( messageCount > MAX_ERRORS ) break;
			}//for
			if( messageCount > MAX_ERRORS) break;
		}//for

		float totalNeighbors = 0.0f;
		for( int particleId = 0; particleId < SPHContstants.PARTICLE_COUNT; ++particleId ){
			float[] myPosition = readField(4 * particleId, solver.sortedPositionPtr );

			for( int neighborNum = 0; neighborNum < SPHContstants.NEIGHBOR_COUNT; ++neighborNum ){
				float neighborParticleId = solver.neighborMapPtr.get(neighborNum * 2);
				if( (int)neighborParticleId == NO_PARTICLE_ID ) continue;
				totalNeighbors += 1.0f;
			}//for
		}//for
		float meanNeighbors = totalNeighbors / SPHContstants.PARTICLE_COUNT;
		if( meanNeighbors < 0.7f * SPHContstants.NEIGHBOR_COUNT ){
			out.println("Error: unexpectedly low meanNeighbors "+meanNeighbors);
			result = false;
		}else{
			out.println("mean neighbors "+meanNeighbors+" NEIGHBOR_COUNT="+SPHContstants.NEIGHBOR_COUNT);
		}

		// truly randomize sampling, assume data is same from run to run
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int particleId = i;
			float[] myPosition = readField( 4 * particleId, solver.sortedPositionPtr );
			for( int neighborParticleId = 0; neighborParticleId < SPHContstants.PARTICLE_COUNT; ++neighborParticleId ){
				if( neighborParticleId == particleId ) continue;
				float[] neighborPosition = readField( 4 * neighborParticleId, solver.sortedPositionPtr );
				float dx = myPosition[ 0 ] - neighborPosition[ 0 ];
				float dy = myPosition[ 1 ] - neighborPosition[ 1 ];
				float dz = myPosition[ 2 ] - neighborPosition[ 2 ];
				float distance = (float)Math.sqrt( dx * dx + dy * dy + dz * dz );
				if( distance < PhysicsConstants.H - EPSILON ){
					boolean found = false;
					for(int k = 0; k < SPHContstants.NEIGHBOR_COUNT; ++k ){
						float npid = solver.neighborMapPtr.get(particleId * SPHContstants.NEIGHBOR_COUNT * 2 + k * 2 );
						if( (int)npid == neighborParticleId ) found = true;
					}//for
					if( !found ){
						out.println("Warning: particle "+particleId+" at ("
							+ myPosition[ 0 ] + ", "
							+ myPosition[ 1 ] +", "
							+ myPosition[ 2 ]+") did not find neighbor "
							+ neighborParticleId +" at ("
							+ neighborPosition[ 0 ] + ", "
							+ neighborPosition[ 1 ] + ", "
							+ neighborPosition[ 2 ]+ ") distance " 
							+ distance);
						messageCount++;
					}
				}
				if( messageCount > MAX_ERRORS ) break;
			}//for
			if( messageCount > MAX_ERRORS ) break;
		}//for
		return result;
	}
		
	private boolean testEndIndexx() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		solver.gridCellIndex.read(solver.queue, solver.gridCellIndexPtr, true);
		solver.queue.finish();

		// test that all values are increasing except for NO_CELL_ID entries
		int lastParticleIndex = 0;
		for( int i = 0; i < solver.gridCellCount + 1; ++i ){
			int particleIndex = solver.gridCellIndexPtr.get(i);
			if( particleIndex < 0 ) continue;
			if( particleIndex < lastParticleIndex ){
				out.println("Error: particle index out of order "+particleIndex+" at grid cell "+i);
				result = false;
			}
			lastParticleIndex = particleIndex;
			if(result == false)
				messageCount++;
			if( messageCount > MAX_ERRORS ) break;
		}//for

		// compute the expected result and compare
		int lastCellId = -1;
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = solver.particleIndexPtr.get(i * 2);
			int particleId = solver.particleIndexPtr.get(i * 2 + 1);
			if( cellId < 0 || cellId > solver.gridCellCount 
				|| particleId < 0 || particleId > SPHContstants.PARTICLE_COUNT ){
					out.println("Error: impossible particle index values, particle "+i
						+" cell "+cellId+" particleId "+particleId);
					result = false;
					if( ++messageCount > MAX_ERRORS )
						break;
					continue;
			}
			if( cellId > lastCellId ){
				int indexValue = solver.gridCellIndexPtr.get(cellId);
				if( indexValue >= 0 && indexValue != i ){
					out.println("Error: grid cell index mismatch at cellId "+cellId+" = " 
						+solver.gridCellIndexPtr.get(cellId)+" expected "+i);
					result = false;
				}
				lastCellId = cellId;
			}
			if( result == false )
				messageCount++;
			if( messageCount > MAX_ERRORS ) break;
		}//for

		// verify the last+1 grid cell index
		if( solver.gridCellIndexPtr.get(solver.gridCellCount) != SPHContstants.PARTICLE_COUNT ){
			out.println("Error: last+1 grid cell index is wrong " 
				+solver.gridCellIndexPtr.get(solver.gridCellCount - 1)+ " should be "+ SPHContstants.PARTICLE_COUNT);
			result = false;
		}
		return result;
	}
	
	private boolean testEndIndexPostPass() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;

		solver.gridCellIndex.read(solver.queue, solver.gridCellIndexPtr, true);
		solver.queue.finish();
		solver.gridCellIndexFixedUp.read(solver.queue, solver.gridCellIndexFixedUpPtr, true);
		solver.queue.finish();
		// test that all values are increasing except for NO_CELL_ID entries
		int[]  _particlesPerCell = new int[ solver.gridCellCount ];
		int lastParticleIndex = 0;
		for( int i = 0; i < solver.gridCellCount + 1; ++i ){
			int particleIndex = solver.gridCellIndexFixedUpPtr.get(i);
			if( particleIndex < lastParticleIndex ){
				out.println("Error: gridCellIndexFixedUp out of order " + particleIndex + " at grid cell "+i);
				result = false;
			}
			if( i > 0 ){
				_particlesPerCell[ i - 1 ] = particleIndex - lastParticleIndex;
			}
			lastParticleIndex = particleIndex;
			if( result == false )
				messageCount += 1;
			if( messageCount > MAX_ERRORS ) break;
		}//for

		float sumParticlesPerCell = 0.0f;
		for( int i = 0; i < solver.gridCellCount; ++i ){
			sumParticlesPerCell += _particlesPerCell[ i ];
		}//for
		float meanParticlesPerCell = sumParticlesPerCell / (float)solver.gridCellCount;
		float totalVariance = 0.0f;
		for( int i = 0; i < solver.gridCellCount; ++i ){
			float varianceThisCell = ( _particlesPerCell[ i ] - meanParticlesPerCell ) * ( _particlesPerCell[ i ] - meanParticlesPerCell );
			totalVariance += varianceThisCell;
		}//for
		float stdDev = (float)Math.sqrt( totalVariance / (float)solver.gridCellCount );

		out.println("Particles per cell, mean " + meanParticlesPerCell + " standard deviation " + stdDev);

		float expectedParticlesPerCell = (float)SPHContstants.PARTICLE_COUNT / (float)solver.gridCellCount;
		if( meanParticlesPerCell != expectedParticlesPerCell ){
			out.println("Error: wrong mean particles per cell " + meanParticlesPerCell 
				+ " should be " + expectedParticlesPerCell);
			result = false;
		}

		float maxExpectedStandardDeviation = (float)expectedParticlesPerCell * 0.06f;//observed empirically 
		// to do - predict expected std dev based on grid size and particle count
		if( stdDev > maxExpectedStandardDeviation ){
			out.println("Warning: nonequilibrium standard deviation " + stdDev 
				+" of particles per cell, equilbrium value would be "+ maxExpectedStandardDeviation );
		}

		// compute the expected result and compare
		solver.particleIndex.read(solver.queue, solver.particleIndexPtr, true);
		solver.queue.finish();
		int lastCellId = -1;
		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			int cellId = solver.particleIndexPtr.get(i * 2);
			int particleId = solver.particleIndexPtr.get(i * 2 + 1);
			if( cellId > lastCellId ){
				int indexValue = solver.gridCellIndexPtr.get(cellId);
				if( indexValue != i ){
					out.println( "Error: grid cell index mismatch at cellId "+ cellId+" = " 
						+solver.gridCellIndexPtr.get(cellId)+" expected " +i );
					result = false;
				}
				lastCellId = cellId;
			}
			if( result == false )
				messageCount += 1;
			if( messageCount > MAX_ERRORS ) break;
		}//for

		// verify the last+1 grid cell index
		if( solver.gridCellIndexPtr.get(solver.gridCellCount) != SPHContstants.PARTICLE_COUNT ){
			out.println("Error: last+1 grid cell index is wrong "+ solver.gridCellIndexPtr.get(solver.gridCellCount - 1)  
			+ " should be "+ SPHContstants.PARTICLE_COUNT );
			result = false;
		}
		return result;
	}
	
	private boolean testEndIntegrate() {
		// TODO Auto-generated method stub
		boolean result = true;
		int messageCount = 0;
		solver.acceleration.read(solver.queue, solver.accelerationPtr, true);
		solver.queue.finish();
		solver.sortedPosition.read(solver.queue,solver.sortedPositionPtr, true);
		solver.queue.finish();
		solver.sortedVelocity.read(solver.queue, solver.sortedVelocityPtr, true);
		solver.queue.finish();
		solver.position.read(solver.queue,solver.positionPtr, true);
		solver.queue.finish();
		solver.velocity.read(solver.queue, solver.velocityPtr, true);
		solver.queue.finish();


		for( int i = 0; i < SPHContstants.PARTICLE_COUNT; ++i ){
			float[] sortedPosition_ = readField( 4 * i, solver.sortedPositionPtr );
			float[] position_ = readField( 4 * i, solver.positionPtr );
			for( int j = 0; j < 3; ++j ){
				if( position_[ 0 ] < SPHContstants.XMIN_F || position_[ 0 ] > SPHContstants.XMAX_F ||
					position_[ 1 ] < SPHContstants.YMIN_F || position_[ 1 ] > SPHContstants.YMAX_F ||
					position_[ 2 ] < SPHContstants.ZMIN_F || position_[ 2 ] > SPHContstants.ZMAX_F ){
						out.println("Error: particle " + i + " escaped the environment at ("
							+ position_[ 0 ] + "," + position_[ 1 ]
						+ "," + position_[ 2 ] + ")");
						result = false;
						break;
				}
			}//for
			if( result == false )
				messageCount += 1;
			if( messageCount > MAX_ERRORS ) break;
		}//for
		return result;
	}

	private float[] readField(int index, Pointer<Float> _field){
		float[] arr = {_field.get(index),_field.get(index + 1), _field.get(index + 2)};
		return arr;
	}
}
