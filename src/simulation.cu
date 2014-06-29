#include "simulation.h" 

__device__ void launch(NeutronInfoStruct* pInfo,int id, CMPTYPE energy){
  pInfo[id].energy = energy;
}

__global__ void initialize(MemStruct pInfo, CMPTYPE energy){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, 1234, 0, &(pInfo.nInfo[id].rndState));
  launch(pInfo.nInfo, id, energy);
  //pInfo[id].energy = energy; //id+1.0; //(id + 1)*1.63*energy*0.001;// 
  pInfo.thread_active[id] = 1u;
  pInfo.tally[id].cnt = 0;

}
#if defined(__QUICKW)
__global__ void initialize_table(CComplex<CMPTYPE> *table){
  int id = blockDim.x*blockIdx.x + threadIdx.x;
  fill_w_tabulated(table, id);
}
#endif

#if defined(__TRACK)
__global__ void history(multipole U238, CMPTYPE* devicearray, MemStruct Info, unsigned num_src, unsigned devstep){
#else
__global__ void history(multipole U238, MemStruct Info, unsigned num_src, unsigned devstep){
#endif
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned istep;
  unsigned live;
  CMPTYPE localenergy;
  CMPTYPE rnd;
  CMPTYPE sigT, sigA, sigF;
  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.nInfo[id].rndState;

  localenergy = Info.nInfo[id].energy;
  unsigned cnt = 0u;
  unsigned terminated = 0u;
  live = 1u;
  //while(live){
  for (istep = 0; istep < devstep; istep++){
    rnd = curand_uniform(&localState);
#if defined(__TRACK)
    unsigned M = gridDim.x*blockDim.x;
    live = Info.tally[id].cnt + cnt;
    live = live*(live<M) + M*(live>=M); 
    if(0==id)
      devicearray[4*live  ] = localenergy;
    if(1==id)
      devicearray[4*live+1] = localenergy;  
    if(2==id)
      devicearray[4*live+2] = localenergy;
    if(3==id)
      devicearray[4*live+3] = localenergy;
#endif

#if defined(__SAMPLE)
    U238.xs_eval_fast(localenergy + 
		      curand_normal(&localState)*sqrt(300.0*KB)*sqrt(0.5)/U238.dev_doubles[SQRTAWR], 
		      sigT, sigA, sigF);
#else
    U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
#endif
    localenergy = localenergy * rnd;
    live = (localenergy > 1.0);
    cnt = cnt + 1;
    /*So far, energy is the only state*/
    localenergy = localenergy*live + 20000.0*(1u - live);
    terminated += !live;
  }
  //}
  /*Note: from now on, live does not indicate neutron but thread active */
  //live = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  atomicAdd(Info.num_terminated_neutrons,terminated);
  Info.thread_active[id] =  (terminated+1)*blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  Info.nInfo[id].rndState = localState; 
  Info.nInfo[id].energy = localenergy;
  Info.tally[id].cnt += cnt; 

}


__global__ void remaining(multipole U238, CMPTYPE *devicearray, MemStruct Info){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned live = true;
  CMPTYPE localenergy;
  CMPTYPE rnd;
  CMPTYPE sigT, sigA, sigF;
 
  /* Copy state to local memory for efficiency */
  curandState localState = Info.nInfo[id].rndState;
  
#if defined(__PROCESS)
  localenergy = 1.0+19999.0/65536.0*id+0.181317676432466;
#else
  localenergy = Info.nInfo[id].energy;
#endif
  unsigned cnt = 0u;
  unsigned terminated = 0u;
  live = 1u;
  while(live){
    rnd = curand_uniform(&localState);
#if defined(__TRACK)
    unsigned M = gridDim.x*blockDim.x;
    live = Info.tally[id].cnt + cnt;
    live = live*(live<M) + M*(live>=M); 
    if(0==id)
      devicearray[4*live  ] = localenergy;
    if(1==id)
      devicearray[4*live+1] = localenergy;  
    if(2==id)
      devicearray[4*live+2] = localenergy;
    if(3==id)
      devicearray[4*live+3] = localenergy;
#endif
#if defined(__SAMPLE)
    U238.xs_eval_fast(localenergy + 
		      curand_normal(&localState)*sqrt(300.0*KB)*sqrt(0.5)/U238.dev_doubles[SQRTAWR], 
		      sigT, sigA, sigF);
#else
    U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
#endif

#if !defined(__PROCESS)
    localenergy = localenergy * rnd;
    live = (localenergy > 1.0);
    cnt = cnt + 1;
    terminated += !live;
#else
    live = false;
#endif
  }
  /* Copy state back to global memory */
  atomicAdd(Info.num_terminated_neutrons,terminated);
  Info.nInfo[id].rndState = localState;
  Info.tally[id].cnt += cnt;

#if !defined(__TRACK)
#if defined(__PROCESS)  
  devicearray[4 * id] = localenergy ;
#else
  devicearray[4 * id] = localenergy / rnd;
#endif
  devicearray[4 * id + 1] = sigT;
  devicearray[4 * id + 2] = sigA;
  devicearray[4 * id + 3] = sigF;
#endif
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

