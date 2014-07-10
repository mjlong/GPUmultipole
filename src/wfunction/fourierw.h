#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"

#define taom 12.0
#define M    23
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z);


#endif
