#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"

#if defined(__CFLOAT)
/* taom and N
  12.0     23
   9.0     13
  8.31     11
  7.52      9
*/
#define taom 8.31f
#else
#define taom 8.31
#endif

#define M    11
__global__ void fill_a(CMPTYPE* a, CMPTYPE* b);
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z);


#endif
