#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

#define STARTENE 20000.0

//TODO: it has not been determined how to save neutron and simulation state

typedef struct {
  CMPTYPE energy;
}basicneutronInfo;

typedef struct {
  unsigned cnt;
  //CMPTYPE   *unknown;
}TallyStruct;

typedef struct {
  CMPTYPE energy;
  curandState rndState;
}NeutronInfoStruct;

typedef struct {
  CMPTYPE energy;
  CMPTYPE sigT;
  CMPTYPE sigA;
  CMPTYPE sigF;
}XsStruct;

typedef struct {
  NeutronInfoStruct *nInfo;
  unsigned int *thread_active;
  unsigned int *num_terminated_neutrons;
  XsStruct *sigma;
  TallyStruct *tally;
}MemStruct;

#if defined (__QUICKW)
#include "QuickW.h"
__global__ void initialize_table(CComplex<CMPTYPE>*);
#endif

#if defined(__TRACK)
__global__ void history(multipole, CMPTYPE *, MemStruct, unsigned, unsigned);
#else
__global__ void history(multipole, MemStruct, unsigned, unsigned );
#endif
__global__ void remaining(multipole, CMPTYPE *, MemStruct );
__global__ void initialize(MemStruct, CMPTYPE);
__device__ void launch(NeutronInfoStruct*, int, CMPTYPE);
__global__ void statistics(TallyStruct*, unsigned*);

#endif
