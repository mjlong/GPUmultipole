#ifndef __FOURIER_H__
#define __FOURIER_H__


#include <stdio.h>
#include "global.h"
#include "CComplex.h"

__device__ CComplex<CMPTYPE> fourierw(CComplex<CMPTYPE> z, CMPTYPE taom, unsigned M);


#endif
