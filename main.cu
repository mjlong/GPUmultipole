#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

__global__ void history(curandState *rndStates, unsigned *);



void anyvalue(struct multipoledata data, int *value, double *d1, double *d2){
  curandState *rndStates;
  unsigned gridx, gridy, blockx, blocky, blockz, gridsize;
  unsigned *hostarray, *devicearray;
  gridx = 3;
  gridy = 2;
  blockx = 4;
  blocky = 5;
  blockz = 6;
  dim3 dimBlock(gridx, gridy);
  dim3 dimGrid(blockx, blocky, blockz);
  gridsize = gridx*gridy*blockx*blocky*blockz;
  cudaMalloc((void**)&rndStates, gridsize*sizeof(curandState));
  cudaMalloc((void**)&devicearray, gridsize*sizeof(unsigned int));
  hostarray = (unsigned*)malloc(gridsize*sizeof(unsigned int));
  multipole U238(data); //host multipoledata to device
  history<<<dimBlock, dimGrid>>>(rndStates, devicearray);
  cudaMemcpy(hostarray, devicearray, gridsize*sizeof(unsigned int), cudaMemcpyDeviceToHost);
  for(int i=0;i<gridsize;i++)
    printf("%5d  %5d\n",i,hostarray[i]);
  return;
}


__global__ void history(curandState *rndStates, unsigned *devicearray){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = 
    (blockDim.x*blockDim.y*blockDim.z)*
    (blockIdx.y*gridDim.x+blockIdx.x)+
    (blockDim.x*blockDim.y)*threadIdx.z+
    blockDim.x*threadIdx.y+
    threadIdx.x;
  bool live=true;
  double energy = 1.0;
  double rnd = 0.0;
  while(live){
    energy = energy * rnd;
    live = (energy>1.0e-4);
  }
  devicearray[id] = id;
}
