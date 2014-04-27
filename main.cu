#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"

__global__ void history();



void anyvalue(struct multipoledata data, int *value, double *d1, double *d2){
  multipole U238(data);
  CPUComplex z;
  cudaMemcpy(value, U238.l_value+25, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&z  , U238.mpdata+10, 2*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(d2, U238.pseudo_rho+1, sizeof(double), cudaMemcpyDeviceToHost);
  *d1 = imag(z);

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
