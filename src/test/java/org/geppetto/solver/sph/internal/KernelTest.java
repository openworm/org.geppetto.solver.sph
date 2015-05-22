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

import static java.lang.System.out;
import java.io.IOException;

import org.bridj.Pointer;
import org.geppetto.solver.sph.SPHSolverService;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem;
import com.nativelibs4java.opencl.CLPlatform.DeviceFeature;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.util.IOUtils;

public class KernelTest {
	private static final float[] TEST_DATA = { 0f, 1f, 2f, 3f, 4f, 5f, 6f, 7f, 8f, 9f };

	/*
	 * TEST JavaCL with CPU
	 * */
	@Test
	public void testKernelCPU() throws Exception {
		out.println("Testing CPU using host memory");
		test(DeviceFeature.CPU, true);
		out.println("Testing CPU using device memory");
		test(DeviceFeature.CPU, false);
	}

	/*
	 * TEST JavaCL with GPU
	 * */
	@Test
	public void testKernelGPU() throws Exception {
		out.println("Testing GPU using host memory");
		test(DeviceFeature.GPU, true);
		out.println("Testing GPU using device memory");
		test(DeviceFeature.GPU, false);
	}

	private void test(DeviceFeature device, boolean useHostMemory) throws Exception {
		FakeSolver fake = new FakeSolver(TEST_DATA, device);

		float[] results;
		if (useHostMemory) {
			results = fake.solveWithHostMemory();
		} else {
			results = fake.solveWithDeviceMemory();
		}

		assertEquals(TEST_DATA.length, results.length, 0);

		for (int i = 0; i < TEST_DATA.length; i++) {
			// Check that every element in the input buffer is equal to
			// the same element in the output buffer
			assertEquals(TEST_DATA[i], results[i], 0);
		}
	}

	public class FakeSolver {
		private CLContext _context;
		private CLQueue _queue;
		private CLKernel _testKernel;
		private float[] _input;
		
		public FakeSolver(float[] input, DeviceFeature device) throws IOException {
			this._input = input;
			
			initializeCL(device);
		}

		private void initializeCL(DeviceFeature feature) throws IOException {
			_context = JavaCL.createBestContext(feature);

			out.println("Created OpenCL context" + _context);
			// an array with available devices
			CLDevice[] devices = _context.getDevices();

			for (int i = 0; i < devices.length; i++) {
				out.println("Found device - " + i + ": " + devices[i]);
			}

			// have a look at the output and select a device
			CLDevice device = devices[0];
			out.println("Using device: " + device);
			out.println("OpenCL version: " + device.getOpenCLVersion());
			out.println("Driver version: " + device.getDriverVersion());
			out.println("Max workgroup size: " + device.getMaxWorkGroupSize());
			out.println("Max workitems size: " + device.getMaxWorkItemSizes()[0]);

			// create command queue on selected device.
			_queue = _context.createDefaultQueue();

			// load sources, create and build program
			String src = IOUtils.readText(SPHSolverService.class.getResourceAsStream("/testKernel.cl"));
			CLProgram program = _context.createProgram(src);

			_testKernel = program.createKernel("test_kernel");
		}

		public synchronized float[] solveWithHostMemory() {
			// Copy input data to a new allocation of native (host) memory
			Pointer<Float> ptrIn = Pointer.pointerToFloats(_input).order(_context.getByteOrder());
			// Allocate native (host) memory for the output data
			Pointer<Float> ptrOut = Pointer.allocateFloats(_input.length).order(_context.getByteOrder());

			// Create CLBuffer wrappers around input/output buffers, no copying
			CLBuffer<Float> bufIn = _context.createBuffer(CLMem.Usage.Input, ptrIn, false);
			CLBuffer<Float> bufOut = _context.createBuffer(CLMem.Usage.Output, ptrOut, false);

			// Setup the method arguments for the kernel
			_testKernel.setArg(0, bufIn);
			_testKernel.setArg(1, bufOut);
			_testKernel.setArg(2, _input.length);

			// Enqueue execution of the kernel
			CLEvent completion = _testKernel.enqueueNDRange(_queue, new int[] { _input.length });
			// Wait for the kernel to execute
			completion.waitFor();

			// Map the output CLBuffer so we can safely read from it
			bufOut.map(_queue, CLMem.MapFlags.Read);
			// Copy output native (host) memory to a JVM-managed float array. This
			// implicitly copies device memory to host memory first, if needed
			float[] output = ptrOut.getFloats();
			// Unmap the output CLBuffer now that reading is finished
			bufOut.unmap(_queue, ptrOut);

			return output;
		}

		public synchronized float[] solveWithDeviceMemory() {
			// Allocate native (device) memory for the input data
			CLBuffer<Float> bufIn = _context.createFloatBuffer(CLMem.Usage.Input, _input.length);
			// Allocate native (device) memory for the output data
			CLBuffer<Float> bufOut = _context.createFloatBuffer(CLMem.Usage.Output, _input.length);

			// Copy input data directly to device memory
			Pointer<Float> ptrIn = bufIn.map(_queue, CLMem.MapFlags.Write);
			ptrIn.setFloats(_input);
			bufIn.unmap(_queue, ptrIn);

			// Setup the method arguments for the kernel
			_testKernel.setArg(0, bufIn);
			_testKernel.setArg(1, bufOut);
			_testKernel.setArg(2, _input.length);

			// Enqueue execution of the kernel
			CLEvent completion = _testKernel.enqueueNDRange(_queue, new int[] { _input.length });

			// Do an explicit copy of the device output buffer to native (host) memory
			// when the kernel finishes executing
			Pointer<Float> output = bufOut.read(_queue, completion);
			// Copy native (host) memory to a JVM-managed float array
			return output.getFloats();
		}
	}
}
