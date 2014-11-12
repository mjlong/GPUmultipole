#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"
#include "gpuerrchk.h"

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

extern __constant__ CMPTYPE a[M+1];
extern __constant__ CMPTYPE b[M+1];

__global__ void fill_a(CMPTYPE* a, CMPTYPE* b);
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z, CMPTYPE* a, CMPTYPE* b);
void release_wtables(CMPTYPE* da, CMPTYPE* db);

__global__ void fill_a(CMPTYPE *a, CMPTYPE *b){
  // since n is supposed to be 1,2,3,...,23
  // by default dimGrid=(1,1,1)
  int n = threadIdx.x+1;
  a[n] = exp(-n*n*PI*PI/(taom*taom));
  b[n] = (CMPTYPE)n*n*PI*PI;
}

__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
  CComplex<CMPTYPE> A = taom*z;//CComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0);//
  CComplex<CMPTYPE> B = exp(ONEI*A);      //CComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0);//
  CComplex<CMPTYPE> C = A*A;
  CComplex<CMPTYPE> w;
  w = CComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0); 
  for(int n=1;n<M;n++){
//    //w = w + exp(-n*n*PI*PI/(taom*taom))*( (CMPTYPE)(((n<<31)>>31)|0x00000001)*B - (CMPTYPE)1.0)/(n*n*PI*PI - C);
//    w = w + a[n]*( -B - (CMPTYPE)1.0)/(b[n] - C);
    n++;
//    w = w + a[n]*(  B - (CMPTYPE)1.0)/(b[n] - C);
  }
////  w = w + a[M] *( -B - (CMPTYPE)1.0)/(b[M]- C);
//
//  w = (w+w)*A;
//  w = w + ((CMPTYPE)1.0-B)/A;
  return w;//*ONEI;
}

#endif
