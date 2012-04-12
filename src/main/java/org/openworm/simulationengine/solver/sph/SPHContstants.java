package org.openworm.simulationengine.solver.sph;

public class SPHContstants {


	public static final int PARTICLE_COUNT = ( 32 * 1024 );
	public static final int NEIGHBOR_COUNT = 32;
	public static final int NK = NEIGHBOR_COUNT * PARTICLE_COUNT;
	
	public static final int XMIN = 0;
	public static final int XMAX = 100;
	public static final int YMIN = 0;
	public static final int YMAX = 40;
	public static final int ZMIN = 0;
	public static final int ZMAX = 40;

	public static final float XMIN_F = XMIN;
	public static final float XMAX_F = XMAX;
	public static final float YMIN_F = YMIN;
	public static final float YMAX_F = YMAX;
	public static final float ZMIN_F = ZMIN;
	public static final float ZMAX_F = ZMAX;
	
}
