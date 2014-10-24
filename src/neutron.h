#ifndef __NEUTRON_H__
#define __NEUTRON_H__

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
  unsigned *isotope;
  CMPTYPE *energy;
  float* pos_x;
  float* pos_y;
  float* pos_z;
  float* dir_polar;
  float* dir_azim;
  CMPTYPE *isoenergy;
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



#endif
