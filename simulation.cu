#include "simulation.h"

__global__ void initialize(neutronInfo Info, double energy){
  int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  Info.energy[id] = (id + 1)*1.63*energy*0.001;

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

  //size of shared[] is given as 3rd parameter while launching the kernel
  extern __shared__ double shared[];
  
  double *tally = &shared[0];

  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &Info.rndState[id]);

  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.rndState[id];

  bool live=true;
  double energy = Info.energy[id];
  devicearray[7*id]=energy;
  double rnd;
  double sigT, sigA, sigF, sibT, sibA, sibF;
  while(live){
    rnd = curand_uniform(&localState);
    //for test
    U238.xs_eval_fast(energy, sigT, sigA, sigF);
    U238.xs_eval_fast(energy, sqrt(900*KB), sibT, sibA, sibF);
    energy = energy * rnd;
    live = false;
  }
   
  devicearray[7*id+1]=sibT;
  devicearray[7*id+2]=sigT;
  devicearray[7*id+3]=sibA;
  devicearray[7*id+4]=sigA;
  devicearray[7*id+5]=sibF;
  devicearray[7*id+6]=sigF;
  
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


