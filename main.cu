#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"

__global__ void history();



void anyvalue(struct multipoledata data, int *value, double *ptr){
  double *prho, *plvalue, *pwstart, *pwend, *pfit;
  CComplex  *psigTfactor, *pmpdata;
  struct multipoledata *pdata;
  multipole U238(data);
  cudaMemcpy(value, U238.l_value+25, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(ptr  , U238.dev_doubles+SQRTAWR, sizeof(double), cudaMemcpyDeviceToHost);
  //  ptr = U238.pseudo_rho;
  dim3 dimBlock(1);
  dim3 dimGrid(10,10,10);
  history<<<dimBlock, dimGrid>>>();
  return;
}


__global__ void history(){
  bool live=true;
  double energy = 1.0;
  while(live){
    energy = energy * 0.8;
    live = (energy>1.0e-4);
  }
}
