#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "multipole.h"
#include "material.h"
#include <cuda.h>
#include <curand_kernel.h>
#include <cudpp.h>
#include <cudpp_config.h>
#include "neutron.h"

#define STARTENE 20000000.0
#define ENDENERG 0.00001
#define MAXENERGY 30000.0

unsigned search_bin(CMPTYPE energy,float* spectrumbins);
#endif
