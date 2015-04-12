#ifndef __H5_READ_H__
#define __H5_READ_H__
#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"

#if defined(__CFLOAT)
#define CMPTYPE float
#else
#define CMPTYPE double 
#endif

typedef struct{
  CMPTYPE complex[2];
}tuple;

void h5read(char filename[]); 

#endif
