#include "simulation.h"

__global__ void initialize(neutronInfo Info, double energy){
  int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  Info.energy[id] = energy;//(id + 1)*1.63*energy*0.001;

}

__global__ void history(multipole U238, double *devicearray, struct neutronInfo Info){
  int i;

  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;

  int idl = 
    (blockDim.x*blockDim.y)*threadIdx.z+
    blockDim.x*threadIdx.y+
    threadIdx.x;

  int idb = 
    blockIdx.y*gridDim.x+blockIdx.x;
  int blocksize = blockDim.x * blockDim.y * blockDim.z;

  bool live=true;
  double energy;
  double rnd;
  double sigT, sigA, sigF;
  extern __shared__ double shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  double *tally = &shared[0];

  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &Info.rndState[id]);

  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.rndState[id];

  energy = Info.energy[id];
  devicearray[4*id]=energy;
  while(live){
    rnd = curand_uniform(&localState);
    U238.xs_eval_fast(energy, sqrt(300.0*KB), sigT, sigA, sigF);
    energy = energy * rnd;
    live = (energy>1.0);
  }
   
  devicearray[4*id+1]=sigT;
  devicearray[4*id+2]=sigA;
  devicearray[4*id+3]=sigF;
  
  /* Copy state back to global memory */ 
  Info.rndState[id] = localState; 

  /*reduce tally*/
  tally[idl] = (double)(rnd<0.5);
  __syncthreads();
  i = blocksize>>1;
  while(i){
    if(idl<i)
      tally[idl] += tally[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    Info.tally[idb] = tally[0]/blocksize;
  }
}


