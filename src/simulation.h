#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include "material.h"
#include <cuda.h>
#include <curand_kernel.h>
#include <cudpp.h>
#include <cudpp_config.h>
#include "neutron.h"

#define STARTENE 20000000.0
#define ENDENERG 0.00001
#define MAXENERGY 30000.0

__global__ void history(material mat, multipole, MemStruct, unsigned num_src,unsigned active);
__global__ void initialize(MemStruct);
__global__ void reduce_sum_plus(unsigned *threadcnt, unsigned* cnt);
__global__ void reduce_sum_equal(unsigned* thread_active, unsigned* active);
__global__ void update_sort_key(MemStruct DeviceMem, material mat);
__global__ void transport(MemStruct, material, unsigned);
__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id);
__global__ void resurrection(NeutronInfoStruct nInfo, unsigned id);


__global__ void z2w_d(CComplex<CMPTYPE> *pz, CComplex<CMPTYPE> *pw);
#endif
