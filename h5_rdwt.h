#ifndef __H5_READ_H__
#define __H5_READ_H__
#include "CPUComplex.h"
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#include "multipole.h"
using namespace std;



void h5read(multipole & pole, char filename[]); 

#endif
