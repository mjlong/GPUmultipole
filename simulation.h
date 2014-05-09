#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>


//TODO: it has not been determined how to save neutron and simulation state

struct basicneutronInfo{
  double energy;
};

struct neutronInfo{
  double *energy;
  curandState *rndState;
  double *tally;
};

__global__ void history(multipole, double*, struct neutronInfo );
__global__ void initialize(struct neutronInfo, double);


#endif