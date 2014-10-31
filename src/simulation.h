#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include "material.h"
#include <cuda.h>
#include <curand_kernel.h>
#include <cudpp.h>
#include <cudpp_config.h>
#include "neutron.h"

#define STARTENE 20000.0
#define MAXENERGY 30000.0

#if defined(__TRACK)
__global__ void history(int, multipole, CMPTYPE *, MemStruct, unsigned);
#else
__global__ void history(int, multipole, MemStruct, unsigned);
#endif
__global__ void remaining(int, multipole, CMPTYPE *, MemStruct );
__global__ void initialize(MemStruct, CMPTYPE);
__device__ void launch(NeutronInfoStruct, int, CMPTYPE);
__global__ void statistics(unsigned*, unsigned*);
__global__ void transport(MemStruct, material);

#endif
