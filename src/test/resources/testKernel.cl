//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#pragma OPENCL EXTENSION cl_intel_printf : enable

__kernel void test_kernel(__global float * in,
                          __global float * out,
                          int SIZE)
{    
    // get index into global data array
    int id = get_global_id(0);

    // bound check (equivalent to the limit on a 'for' loop for standard/serial C code
    if (id >= SIZE)  {
        return;
    }
    
    out[id] = in[id];
}