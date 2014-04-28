#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

__global__ void history(multipole, curandState *rndState, double *);




void anyvalue(struct multipoledata data, int *value, double *d1, double *d2){
  curandState *rndState;
  unsigned gridx, gridy, blockx, blocky, blockz, gridsize;
  double *hostarray, *devicearray;
  gridx = 40;
  gridy = 40;
  blockx = 5;
  blocky = 5;
  blockz = 5;
  dim3 dimBlock(gridx, gridy);
  dim3 dimGrid(blockx, blocky, blockz);
  gridsize = gridx*gridy*blockx*blocky*blockz;
  cudaMalloc((void**)&rndState, gridsize*sizeof(curandState));
  cudaMalloc((void**)&devicearray, 4*gridsize*sizeof(double));
  hostarray = (double*)malloc(4*gridsize*sizeof(double));
  multipole U238(data); //host multipoledata to device
  history<<<dimBlock, dimGrid>>>(U238, rndState, devicearray);
  cudaMemcpy(hostarray, devicearray, 4*gridsize*sizeof(double), cudaMemcpyDeviceToHost);


  for(int i=0;i<gridsize;i++){
    printf("%9.3e,  %10.6e, %10.6e, %10.6e\n",
	   hostarray[4*i],
	   hostarray[4*i+1],
	   hostarray[4*i+2],
	   hostarray[4*i+3]);
  }

  return;
}


__global__ void history(multipole U238, curandState *rndState, double *devicearray){
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
  double energy = 1000.0;
  double rnd;
  double sigT, sigA, sigF;
  while(live){
    rnd = curand_uniform(&localState);
    energy = energy * rnd;

    sigT = energy*1.2;
    sigA = energy*0.8;
    sigF = energy*0.02;
    U238.xs_eval_fast(energy, sigT, sigA, sigF);

    live = (energy>50.0);
  }

  devicearray[4*id]=energy;//id*1.0;//energy;
  devicearray[4*id+1]=sigT;//U238.dev_doubles[SPACING];
  devicearray[4*id+2]=sigA;//U238.dev_doubles[STARTE];
  devicearray[4*id+3]=sigF;//U238.dev_doubles[SQRTAWR];

  /* Copy state back to global memory */ 
  rndState[id] = localState; 

}
