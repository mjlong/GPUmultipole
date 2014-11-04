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
__global__ void history(material mat, multipole, CMPTYPE *, MemStruct, unsigned);
#else
__global__ void history(material mat, multipole, MemStruct, unsigned);
#endif
__global__ void remaining(material mat, multipole, CMPTYPE *, MemStruct );
__global__ void initialize(MemStruct, CMPTYPE);
__device__ void launch(NeutronInfoStruct, int, CMPTYPE);
__global__ void statistics(unsigned*, unsigned*);
__global__ void update_sort_key(MemStruct DeviceMem, material mat);
__global__ void transport(MemStruct, material);
__device__ void source_sampling(NeutronInfoStruct nInfo, unsigned id);
__global__ void resurrection();
#endif
