#ifndef __FOURIER_H__
#define __FOURIER_H__

#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double 
#define CMPTYPE2 double2
#endif

#include <stdio.h>
#include "global.h"
#include "CComplex.h"

__device__ CComplex<CMPTYPE> fourierw(CComplex<CMPTYPE> z, CMPTYPE taom, unsigned M);


#endif
