#ifndef __PROCESS_H__
#define __PROCESS_H__
#include <stdio.h>
#include "neutron.h"
void getASE(float *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, float ref, float* ASE);
void cnt2flux(MemStruct HostMem, unsigned numhis, float dx, unsigned meshes, unsigned nbat);
float autok(float *batmeans, unsigned n, unsigned k, unsigned meshes, unsigned im);
void getCOR(float *batmeans, unsigned meshes, unsigned nbat, unsigned ubat, unsigned upto,float *COR);
#endif
