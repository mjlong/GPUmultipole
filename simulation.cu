#include "simulation.h" 

__device__ void launch(NeutronInfoStruct* pInfo,int id, double energy){
  pInfo[id].energy = energy;
}

__global__ void initialize(MemStruct pInfo, double energy){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  launch(pInfo.nInfo, id, energy);
  //pInfo[id].energy = energy; //id+1.0; //(id + 1)*1.63*energy*0.001;// 
  pInfo.thread_active[id] = 1u;
  pInfo.tally[id].cnt = 0;

}

__global__ void history(multipole U238, MemStruct Info, unsigned num_src, unsigned devstep){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned istep;
  bool live=true;
  double localenergy,initenergy;
  double rnd;
  double sigT, sigA, sigF;

  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &(Info.nInfo[id].rndState));

  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.nInfo[id].rndState;

  initenergy = Info.nInfo[id].energy;
  localenergy = initenergy;
  unsigned cnt = 0u;
  unsigned terminated = 0u;
  //while(live){
  for (istep = 0; istep < devstep; istep++){
	  rnd = curand_uniform(&localState);
	  U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
	  localenergy = localenergy * rnd;
	  live = (localenergy > 1.0);
	  cnt = cnt + 1;
          if(!live){
            terminated = atomicAdd(Info.num_terminated_neutrons, 1u);
	    localenergy = 2000.0;
            if(terminated >= (num_src - blockDim.x)){
              istep = devstep;  
              Info.thread_active[id] = 0u;
            }
          }
  }
  //}
  Info.nInfo[id].energy = localenergy;
  /* Copy state back to global memory */ 
  Info.nInfo[id].rndState = localState; 
  Info.tally[id].cnt += cnt; 

}

__global__ void remaining(multipole U238, double *devicearray, MemStruct Info){
	//TODO:this is one scheme to match threads to 1D array, 
	//try others when real simulation structure becomes clear
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	bool live = true;
	double localenergy;
	double rnd;
	double sigT, sigA, sigF;

	/* Each thread gets same seed, a different sequence number, no offset */
	curand_init(1234, id, 0, &(Info.nInfo[id].rndState));

	/* Copy state to local memory for efficiency */
	curandState localState = Info.nInfo[id].rndState;

	localenergy = Info.nInfo[id].energy;
	unsigned cnt = 0u;
	while(live){
		rnd = curand_uniform(&localState);
		U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
		localenergy = localenergy * rnd;
		live = (localenergy > 1.0);
		cnt = cnt + 1;
	}
	/* Copy state back to global memory */
	Info.nInfo[id].rndState = localState;
	Info.tally[id].cnt += cnt;

	devicearray[4 * id] = localenergy / rnd;
	devicearray[4 * id + 1] = sigT;
	devicearray[4 * id + 2] = sigA;
	devicearray[4 * id + 3] = sigF;
}


__global__ void statistics(TallyStruct *threadtally, unsigned* cnt){
  /*reduce tally*/
  /*TODO:
    alternatives:
    1. only count for a block, saving global memory (acceess)
    2. count for each thread, saving time in thread wait
  */
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ unsigned shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  shared[idl] = threadtally[id].cnt;
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      shared[idl] += shared[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    cnt[blockIdx.x] = shared[0];
  }
  
}


__global__ void isActive(MemStruct DevMem, unsigned int *active){
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ unsigned shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  shared[idl] = DevMem.thread_active[id]; 
  __syncthreads();
  i = blockDim.x>>1;
  while(i){
    if(idl<i)
      shared[idl] += shared[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    active[blockIdx.x] = shared[0];
  }


}

