#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

__global__ void history(curandState *rndState, double *);



void anyvalue(struct multipoledata data, int *value, double *d1, double *d2){
  curandState *rndState;
  unsigned gridx, gridy, blockx, blocky, blockz, gridsize;
  double *hostarray, *devicearray;
  gridx = 3;
  gridy = 2;
  blockx = 4;
  blocky = 5;
  blockz = 6;
  dim3 dimBlock(gridx, gridy);
  dim3 dimGrid(blockx, blocky, blockz);
  gridsize = gridx*gridy*blockx*blocky*blockz;
  cudaMalloc((void**)&rndState, gridsize*sizeof(curandState));
  cudaMalloc((void**)&devicearray, gridsize*sizeof(double));
  hostarray = (double*)malloc(gridsize*sizeof(double));
  multipole U238(data); //host multipoledata to device
  history<<<dimBlock, dimGrid>>>(rndState, devicearray);
  cudaMemcpy(hostarray, devicearray, gridsize*sizeof(double), cudaMemcpyDeviceToHost);
  double x, s1=0.0, s2=0.0;
  /*
  for(int i=0;i<gridsize;i++){
    x = hostarray[i];
    printf("%5d  %10.6e\n",i,x);
    s1 += x;
    s2 += x*x;
  }
  printf("mean=%3.1f, var*12=%3.1f, var*12=%3.1f\n", 
	 s1/gridsize, 
	 (s2-s1/gridsize*s1)/(gridsize-1.0)*12.0,
	 (s2/gridsize - s1*s1/gridsize/gridsize)*12.0);
  */
  return;
}


__global__ void history(curandState *rndState, double *devicearray){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = 
    (blockDim.x*blockDim.y*blockDim.z)*
    (blockIdx.y*gridDim.x+blockIdx.x)+
    (blockDim.x*blockDim.y)*threadIdx.z+
    blockDim.x*threadIdx.y+
    threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &rndState[id]);
  /* Copy state to local memory for efficiency */ 
  curandState localState = rndState[id];
  bool live=true;
  double energy = 1.0;
  double rnd;
  while(live){
    rnd = curand_uniform(&localState);
    energy = energy * rnd;
    live = (energy>1.0e-4);
  }
  /* Copy state back to global memory */ 
  rndState[id] = localState; 

}
