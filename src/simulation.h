#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

#define STARTENE 20000.0
#define MAXENERGY 30000.0

//TODO: it has not been determined how to save neutron and simulation state
//TODO: array of struct Vs struct of array
typedef struct {
  CMPTYPE energy;
}basicneutronInfo;

typedef struct {
  unsigned *cnt;
  //CMPTYPE   *unknown;
}TallyStruct;

typedef struct {
  unsigned *isotope;
  CMPTYPE *energy;
  curandState *rndState;
}NeutronInfoStruct;

typedef struct {
  CMPTYPE energy;
  CMPTYPE sigT;
  CMPTYPE sigA;
  CMPTYPE sigF;
}XsStruct;

typedef struct {
  NeutronInfoStruct nInfo;
  unsigned int *block_terminated_neutrons;
  unsigned int *num_terminated_neutrons;
  XsStruct *sigma;
  TallyStruct tally;
}MemStruct;


#if defined(__TRACK)
__global__ void history(int, multipole, CMPTYPE *, MemStruct, unsigned, unsigned);
#else
__global__ void history(int, multipole, MemStruct, unsigned, unsigned );
#endif
__global__ void remaining(int, multipole, CMPTYPE *, MemStruct );
__global__ void initialize(MemStruct, CMPTYPE);
__device__ void launch(NeutronInfoStruct, int, CMPTYPE);
__global__ void statistics(unsigned*, unsigned*);

#endif
