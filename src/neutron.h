#ifndef __NEUTRON_H__
#define __NEUTRON_H__

#include <curand_kernel.h>

#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double 
#define CMPTYPE2 double2
#endif

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
  unsigned *id;
  unsigned *imat;
  CMPTYPE *energy;
  CMPTYPE *sigT;
  CMPTYPE *sigA;
  CMPTYPE *sigF;
  float* pos_x;
  float* pos_y;
  float* pos_z;
  float* dir_polar;
  float* dir_azimu;
  float* d_closest;
  CMPTYPE *isoenergy;
  curandState *rndState;
}NeutronInfoStruct;

typedef struct {
  NeutronInfoStruct nInfo;
  unsigned int *block_terminated_neutrons;
  unsigned int *num_terminated_neutrons;
  TallyStruct tally;
}MemStruct;



#endif
