#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"

#if defined(__INTERPEXP)
#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double 
#define CMPTYPE2 double2
#endif

#define LENGTH 62
#define WIDTH  62
#endif
/* taom and N
  12.0     23
   9.0     13
  8.31     11
  7.93     10
  7.52      9
*/
#if defined(__CFLOAT)
#define taom 7.93f
#else
#define taom 7.93
#endif

#define M    10
#if defined(__INTERPEXP)
__global__ void fill_exp_table(CComplex<CMPTYPE> *);
#endif
__global__ void fill_a(CMPTYPE* a, CMPTYPE* b);
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z);
#if defined(__PFOURIERW)
__device__ CComplex<CMPTYPE> w_part(CComplex<CMPTYPE> z, unsigned n, CMPTYPE coef);
#endif

#endif
