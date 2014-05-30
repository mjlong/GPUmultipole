#include "simulation.h" 
__global__ void initialize(neutronInfo Info, double energy){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  Info.energy[id] = energy; //(id+1.0); //energy;//(id + 1)*1.63*energy*0.001;// 

}

__global__ void history(multipole U238, double *devicearray, struct neutronInfo Info){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;

  bool live=true;
  double localenergy;
  double rnd;
  double sigT, sigA, sigF;
  extern __shared__ unsigned shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &Info.rndState[id]);

  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.rndState[id];

  localenergy = Info.energy[id];
  unsigned cnt = 0;
  while(live){
    rnd = curand_uniform(&localState);
    U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
    localenergy = localenergy * rnd;
    live = (localenergy>1.0);
    cnt = cnt + 1;
    //live = false;
  }
   
  devicearray[4*id]=localenergy/rnd;
  devicearray[4*id+1]=sigT;
  devicearray[4*id+2]=sigA;
  devicearray[4*id+3]=sigF;
  
  /* Copy state back to global memory */ 
  Info.rndState[id] = localState; 

  /*reduce tally*/
  unsigned *tally = shared;
  int idl = threadIdx.x;
  int i;
  tally[idl] = cnt;
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      tally[idl] += tally[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    Info.tally[blockIdx.x] = tally[0];
  }
}


