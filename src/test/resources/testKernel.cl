
//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#pragma OPENCL EXTENSION cl_intel_printf : enable


__kernel void test_kernel(__global float4 * in,
                          __global float4 * out){
    int id = get_global_id(0);
    out[id] = in[id];
}