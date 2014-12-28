#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "neutron.h"

#define STARTENE 20000000.0
#define ENDENERG 0.00001
#define MAXENERGY 30000.0

unsigned search_bin(CMPTYPE energy,float* spectrumbins);

#if defined(__W__GPU)||defined(__XS_GPU)
#include <cuda.h>
#include <curand_kernel.h>
#include <cudpp.h>
#include <cudpp_config.h>
#include "CComplex.h"
#include "Faddeeva.h"
#endif
#if defined(__W__GPU)
#if defined(__MITW)
#define w_function Faddeeva::w
#endif
#if defined(__QUICKW)
#include "QuickW.h"
#endif
#if defined(__FOURIERW)
#include "fourierw.h"
#endif
#include <stdio.h>
__global__ void device_w_eval(CComplex<CMPTYPE>* z_d, CComplex<CMPTYPE>* w_d);
#endif


#if defined(__XS_GPU)
#include "multipole.h"
#include "material.h"

__global__ void device_xs_eval(multipole mp_para, unsigned *iS, CMPTYPE E, CMPTYPE sqrtKT, 
                               CMPTYPE* sigTs, CMPTYPE* sigAs, CMPTYPE* sigFs);
#endif
#endif
