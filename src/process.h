#ifndef __PROCESS_H__
#define __PROCESS_H__
#include <stdio.h>
#include "neutron.h"
void getASE(double *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, double ref, double* ASE);
void cnt2flux(MemStruct HostMem, unsigned numhis, float dx, unsigned meshes, unsigned nbat);
double autok(double *batmeans, unsigned n, unsigned k, unsigned meshes, unsigned im);
void getCOR(double *batmeans, unsigned meshes, unsigned nbat, unsigned ubat, unsigned upto,double *COR);
#endif
