#ifndef __MATERIALDATA_H__
#define __MATERIALDATA_H__

#include <stdio.h>
#include <stdlib.h>

#define MAXLEN 128
struct matdata{
  unsigned numMat;
  unsigned* offsets;
  float* N_tot; //total number density
  float* densities;//number density
  unsigned *isotopes;
};

void clearline(char *line);
unsigned matread(struct matdata*, char*);
#endif
