#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"

#if defined(__CFLOAT)
#define taom 9.0f //12.0f
#else
#define taom 9.0  //12.0
#endif

#define M    13   //23
__global__ void fill_a(CMPTYPE* a, CMPTYPE* b);
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z);


#endif
