#ifndef __H5_READ_H__
#define __H5_READ_H__
#include <complex>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#include "multipole.h"
using namespace std;

typedef struct{
  double complex[2];
}tuple;


void h5read(multipole & pole, char filename[]); 

#endif
