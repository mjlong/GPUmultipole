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

#define NUM_BINS 16
//TODO: it has not been determined how to save neutron and simulation state
//TODO: array of struct Vs struct of array
typedef struct {
  CMPTYPE energy;
}basicneutronInfo;

typedef struct {
  unsigned cnt;
  //CMPTYPE   *unknown;
}TallyStruct;

typedef struct {
  unsigned imat;
  unsigned live;
  CMPTYPE energy;
  CMPTYPE sigT;
  CMPTYPE sigA;
  CMPTYPE sigF;
}NeutronInfoStruct;

typedef struct {
  NeutronInfoStruct nInfo;
  unsigned num_terminated_neutrons;
  unsigned int *spectrum;
  float *tallybins;
  TallyStruct tally;
}MemStruct;



#endif
