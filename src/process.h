#ifndef __PROCESS_H__
#define __PROCESS_H__
#include "neutron.h"
void getASE(float *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, float ref, float* ASE);
void cnt2flux(MemStruct HostMem, unsigned numhis, float dx, unsigned meshes, unsigned nbat);
#endif
