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

//fill char line[] with ' ' to clear history
void clearline(char *line);

//read materials (isotope and number density) from file
unsigned matread(struct matdata*, char* file);

//after copy data froom matdata to device, release memory
void freeMaterialData(struct matdata* pdata);
#endif
