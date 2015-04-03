#ifndef __PROCESS_H__
#define __PROCESS_H__
#include <stdio.h>
#include "neutron.h"
void getASE(double *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, double ref, double* ASE);
void cnt2flux(MemStruct HostMem, unsigned numhis, float dx, unsigned meshes, unsigned nbat);
double autok(double *batmeans, unsigned n, unsigned k, unsigned meshes, unsigned im);
void getCOR(double *batmeans, unsigned meshes, unsigned nbat, unsigned ubat, unsigned upto,double *COR);
void fitall(double *rhos,unsigned upto, unsigned meshes, double *rho0s, double *qs);
void fitall1(double *rhos,unsigned upto, unsigned meshes, double *rho0s, double *qs);
void fitrho1(double* rho, unsigned m, double* rho0, double* q);
void fitrho(double* rho, unsigned m, double* rho0, double* q);
double variance(double *x, unsigned n);
double variance(double *x, unsigned nbat,unsigned ubat, unsigned meshes, unsigned im);
void varcorrect(double rho0,double q,unsigned m, double *correct);
#endif
