package org.openworm.simulationengine.solver.sph.internal;

import static java.lang.System.out;
import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.bridj.Pointer;
import org.junit.After;
import org.junit.Test;
import org.openworm.simulationengine.solver.sph.SPHSolverService;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.opencl.CLPlatform.DeviceFeature;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.util.IOUtils;

public class KernelTest
{

	private CLContext _context;
	public CLQueue _queue;
	private CLProgram _program;
	private CLDevice _device;
	public CLBuffer<Float> _buffer1;
	public CLBuffer<Float> _buffer2;
	private CLKernel _test_kernel;
	private Pointer<Float> _buffer1Ptr;
	private Pointer<Float> _buffer2Ptr;
	private static final int SIZE=10;
	
	private void init(DeviceFeature device) throws IOException  
	{
		_context = JavaCL.createBestContext(device);
		
		out.println("created "+ _context);
		// an array with available devices
		CLDevice[] devices = _context.getDevices();

		for(int i=0; i<devices.length; i++)
		{
			out.println("device - " + i + ": " + devices[i]);
		}	

		// have a look at the output and select a device
		_device = devices[0];
		out.println("Version " + _device.getOpenCLVersion());
		out.println("Version " + _device.getDriverVersion());
		out.println("using "+ _device);
		out.println("max workgroup size: " + _device.getMaxWorkGroupSize());
		out.println("max workitems size: " + _device.getMaxWorkItemSizes()[0]);
		
		// create command queue on selected device.
		_queue = _context.createDefaultQueue();
		
		// load sources, create and build program
		String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/testKernel.cl"));
		_program = _context.createProgram(src);
		
	 	_test_kernel = _program.createKernel("test_kernel");
	 	
		_buffer1Ptr = Pointer.allocateFloats(SIZE);
		_buffer2Ptr = Pointer.allocateFloats(SIZE);
	 	_buffer1 = _context.createBuffer(Usage.InputOutput,_buffer1Ptr,false);
	 	_buffer2 = _context.createBuffer(Usage.InputOutput, _buffer2Ptr,false);
	 	
	 	_queue.finish();
	}

	public int runTestKernel(){
	 	// Stage SortPostPass
	 	_test_kernel.setArg( 0, _buffer1 );
	 	_test_kernel.setArg(1, _buffer2 );
	 	_test_kernel.enqueueNDRange(_queue, new int[] {SIZE});
	 	return 0;
	 }
	
	private float b1[]={0f,1f,2f,3f,4f,5f,6f,7f,8f,9f};
	private float b2[]={0f,0f,0f,0f,0f,0f,0f,0f,0f,0f};
	
	private void test(DeviceFeature device) throws Exception
	{
		init(device);
		for(int i=0;i<10;i++)
		{
			_buffer1Ptr.set(i,b1[i]);
			_buffer2Ptr.set(i,b2[i]);
		}
		
		//the kernel copies the first buffer in the second one
		runTestKernel();
		
		_buffer1Ptr = _buffer1.read(_queue);
		_buffer2Ptr = _buffer2.read(_queue);
		_queue.finish();
		
		for(int i=0;i<10;i++)
		{
			//so we check that every element in the second buffer is equal to the same element in the first buffer
			assertEquals(_buffer2Ptr.get(i).floatValue(),b1[i],0);
		}
		
		
	}
	
	@Test
	public void testKernelCPU() throws Exception
	{
		test(DeviceFeature.CPU);
	}
	
	@Test
	public void testKernelGPU() throws Exception
	{
		test(DeviceFeature.GPU);
	}
	
	@After
	public void cleanup()
	{
		System.out.println("Cleanup");
		_queue.finish();
		_context.release();
	}
	
}
