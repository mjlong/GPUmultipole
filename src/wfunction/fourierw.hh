#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CPUComplex.h"

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
void fill_a(CMPTYPE* a, CMPTYPE* b);
CPUComplex<CMPTYPE> w_function(CPUComplex<CMPTYPE> z,CMPTYPE*, CMPTYPE*);


#endif
