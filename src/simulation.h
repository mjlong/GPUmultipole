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
__global__ void history(MemStruct, unsigned num_src,unsigned active,unsigned devstep,float width, float dx);
__global__ void initialize(MemStruct,float);
__global__ void reduce_sum_plus(unsigned *threadcnt, unsigned* cnt);
__global__ void reduce_sum_equal(unsigned* thread_active, unsigned* active);
#endif
