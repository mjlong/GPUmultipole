#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

#define NUMSRC 1000000

//TODO: it has not been determined how to save neutron and simulation state

typedef struct {
  double energy;
}basicneutronInfo;

typedef struct {
  unsigned cnt;
  //double   *unknown;
}TallyStruct;

typedef struct {
  double energy;
  curandState rndState;
}NeutronInfoStruct;

typedef struct {
  double energy;
  double sigT;
  double sigA;
  double sigF;
}XsStruct;

typedef struct {
  NeutronInfoStruct *nInfo;
  unsigned int *thread_active;
  XsStruct *sigma;
  TallyStruct *tally;
}MemStruct;


__global__ void history(multipole, double*, NeutronInfoStruct*, TallyStruct* );
__global__ void initialize(MemStruct, double);
__device__ void launch(NeutronInfoStruct*, int, double);
__global__ void statistics(TallyStruct*, unsigned*);

#endif
