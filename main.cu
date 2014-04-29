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
  gridx = 5;
  gridy = 2;
  blockx = 1;
  blocky = 1;
  blockz = 1;
  dim3 dimBlock(gridx, gridy);
  dim3 dimGrid(blockx, blocky, blockz);
  gridsize = gridx*gridy*blockx*blocky*blockz;
  cudaMalloc((void**)&rndState, gridsize*sizeof(curandState));
  cudaMalloc((void**)&devicearray, 7*gridsize*sizeof(double));
  hostarray = (double*)malloc(7*gridsize*sizeof(double));
  multipole U238(data); //host multipoledata to device
  history<<<dimBlock, dimGrid>>>(U238, rndState, devicearray);
  cudaMemcpy(hostarray, devicearray, 7*gridsize*sizeof(double), cudaMemcpyDeviceToHost);


  for(int i=0;i<gridsize;i++){
    printf("%8.4f %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n",
	   hostarray[7*i],
	   hostarray[7*i+1],
	   hostarray[7*i+2],
	   hostarray[7*i+3],
	   hostarray[7*i+4],
	   hostarray[7*i+5],
	   hostarray[7*i+6]);
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
  double sigT, sigA, sigF, sibT, sibA, sibF;
  while(live){
    rnd = curand_uniform(&localState);
    energy = energy * rnd;
    //for test
    energy = (1+id)*1.63;
    U238.xs_eval_fast(energy, sigT, sigA, sigF);
    U238.xs_eval_fast(energy, sqrt(900*KB), sibT, sibA, sibF);

    live = false;//(energy>50.0);
  }

  devicearray[7*id]=energy;
  devicearray[7*id+1]=sibT;
  devicearray[7*id+2]=sigT;
  devicearray[7*id+3]=sibA;
  devicearray[7*id+4]=sigA;
  devicearray[7*id+5]=sibF;
  devicearray[7*id+6]=sigF;

  /* Copy state back to global memory */ 
  rndState[id] = localState; 

}
