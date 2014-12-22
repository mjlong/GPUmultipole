#ifndef __H5_READ_H__
#define __H5_READ_H__
#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#include "CPUComplex.h"
#include "multipole.h"
//#include "multipole_data.h"

#define CMPTYPE double 

typedef struct{
  CMPTYPE complex[2];
}tuple;

void h5read(multipole & pole, char filename[]); 

#endif

