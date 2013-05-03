package org.geppetto.solver.sph.internal;

import static java.lang.System.out;
import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.ByteOrder;

import org.bridj.Pointer;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.After;
import org.junit.Test;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem.MapFlags;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.opencl.CLPlatform.DeviceFeature;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.util.IOUtils;

public class KernelTest {
	
	private int buff_size = 10;
	private float buff1[] = { 0f, 1f, 2f, 3f, 4f, 5f, 6f, 7f, 8f, 9f };
	private float buff2[] = { 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f };

	private void test(DeviceFeature device) throws Exception {
		FakeSolver fake = new FakeSolver(buff1, buff2, buff_size, device);
		
		float[] buff2_results = fake.solve();
		
		for (int i = 0; i < buff_size; i++) {
			// so we check that every element in the second buffer is equal to
			// the same element in the first buffer
			assertEquals(buff1[i], buff2_results[i], 0);
		}
	}

	/*
	 * TEST JavaCL with CPU
	 * */
	@Test
	public void testKernelCPU() throws Exception {
		test(DeviceFeature.CPU);
	}
	
	/*
	 * TEST JavaCL with GPU
	 * */
	@Test
	public void testKernelGPU() throws Exception {
		//test(DeviceFeature.GPU);
	}

	public class FakeSolver{
		private CLContext _context;
		public CLQueue _queue;
		private CLProgram _program;
		private CLDevice _device;
		private CLKernel _test_kernel;
		public CLBuffer<Float> _buffer1;
		public CLBuffer<Float> _buffer2;
		private Pointer<Float> _buffer1Ptr;
		private Pointer<Float> _buffer2Ptr;
		private int SIZE;
		
		private float b1[] = null;
		private float b2[] = null;
		
		public FakeSolver(float[] buff1, float[] buff2, int size, DeviceFeature device) throws IOException {
			this.SIZE = size;
			this.b1 = buff1;
			this.b2 = buff2;
			
			initializeCL(device);
		}
		
		private void initializeCL(DeviceFeature device) throws IOException {
			_context = JavaCL.createBestContext(device);
			ByteOrder byteOrder = _context.getByteOrder();

			out.println("created " + _context);
			// an array with available devices
			CLDevice[] devices = _context.getDevices();

			for (int i = 0; i < devices.length; i++) {
				out.println("device - " + i + ": " + devices[i]);
			}

			// have a look at the output and select a device
			_device = devices[0];
			out.println("Version " + _device.getOpenCLVersion());
			out.println("Version " + _device.getDriverVersion());
			out.println("using " + _device);
			out.println("max workgroup size: " + _device.getMaxWorkGroupSize());
			out.println("max workitems size: " + _device.getMaxWorkItemSizes()[0]);

			// create command queue on selected device.
			_queue = _context.createDefaultQueue();

			// load sources, create and build program
			String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/testKernel.cl"));
			_program = _context.createProgram(src);

			_test_kernel = _program.createKernel("test_kernel");

			// example with host pointer allocation
			_buffer1Ptr = Pointer.allocateFloats(SIZE).order(byteOrder);
			_buffer2Ptr = Pointer.allocateFloats(SIZE).order(byteOrder);
			_buffer1 = _context.createBuffer(Usage.Input, _buffer1Ptr, false);
			_buffer2 = _context.createBuffer(Usage.Output, _buffer2Ptr, false);
			
			// example with direct mapping
			//_buffer1 = _context.createFloatBuffer(Usage.Input, SIZE);
			//_buffer2 = _context.createFloatBuffer(Usage.Output, SIZE);
			//_buffer1Ptr = _buffer1.map(_queue, MapFlags.Write);
     		//_buffer2Ptr = _buffer2.map(_queue, MapFlags.ReadWrite);
			
			// populate buffers
			for (int i = 0; i < SIZE; i++) {
				_buffer1Ptr.set(i, b1[i]);
				_buffer2Ptr.set(i, b2[i]);
			}
		}
		
		public float[] solve(){	
			// the kernel copies the first buffer in the second one
			_test_kernel.setArg(0, _buffer1);
			_test_kernel.setArg(1, _buffer2);
			_test_kernel.setArg(2, SIZE);
			CLEvent testEvt = _test_kernel.enqueueNDRange(_queue, new int[] { SIZE });

			// shouldn't need to do this read operation if the buffers are mapped?
			_buffer2Ptr = _buffer2.read(_queue, testEvt);
			
			float[] results = new float[SIZE];
			for (int i = 0; i < SIZE; i++) {
				results[i] = _buffer2Ptr.get(i);
			}
			
			// un-map buffers if they were mapped
			//_buffer1.unmap(_queue, _buffer1Ptr);
			//_buffer2.unmap(_queue, _buffer2Ptr);
			
			return results;
		}	
	}
}