#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <cuda.h>
#include <curand_kernel.h>
#include "neutron.h"
#include "global.h"

#define STARTENE 20000000.0
#define ENDENERG 0.00001
#define MAXENERGY 30000.0
__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id,float width);
__global__ void history(MemStruct, unsigned num_src,unsigned active,unsigned banksize);
__global__ void initialize(MemStruct,float,int);
__global__ void reduce_sum_plus(int *threadcnt, int* cnt);
__global__ void reduce_sum_equal(int* thread_active, int* active);
__device__ unsigned notleak(float x,float a);
__global__ void history_ref(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize);
__global__ void history_3d_ref(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize);

#define TEPSILON 1.0e-5
#endif
