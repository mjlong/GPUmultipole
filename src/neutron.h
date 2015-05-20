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
  CMPTYPE *cnt;
  //CMPTYPE   *unknown;
}TallyStruct;

typedef struct {
  unsigned *id;
  int *imat;
  int *live;
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

  float* d_pos_x;
  float* d_pos_y;
  float* d_pos_z;
  float* d_time;
  int*   d_igen;
  int*   d_nu;
}NeutronInfoStruct;

typedef struct {
  float *wdspp;
  NeutronInfoStruct nInfo;
  int *block_terminated_neutrons;
  int *num_terminated_neutrons;
  int *num_live_neutrons;
  int *initial_delayed;
  CMPTYPE *block_spectrum;
  CMPTYPE *spectrum;
  double *batchmeans;
  double *accmeans;
  CMPTYPE *batcnt;
  CMPTYPE *tallybins;
  TallyStruct tally;
}MemStruct;



#endif
