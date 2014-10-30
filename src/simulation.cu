#include "simulation.h" 

__device__ void launch(NeutronInfoStruct pInfo,int id, CMPTYPE energy){
  pInfo.energy[id] = energy;
}

__global__ void initialize(MemStruct pInfo, CMPTYPE energy){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &(pInfo.nInfo.rndState[id]));
  curandState state = pInfo.nInfo.rndState[id];
  launch(pInfo.nInfo, id, energy);

//TODO: source sampling should take settings dependent on geometry
  pInfo.nInfo.pos_x[id] = 0.5f+curand_uniform(&state);
  pInfo.nInfo.pos_y[id] = 0.5f+curand_uniform(&state);
  pInfo.nInfo.pos_z[id] = 0.5f+curand_uniform(&state);
  pInfo.nInfo.dir_polar[id] = curand_uniform(&state)*PI;
  pInfo.nInfo.dir_azimu[id] = curand_uniform(&state)*PI*2;

  //pInfo[id].energy = energy; //id+1.0; //(id + 1)*1.63*energy*0.001;// 
  pInfo.nInfo.id[id] = id;
  pInfo.nInfo.isotope[id]=id%2;//0;//
  pInfo.nInfo.isoenergy[id]=MAXENERGY*(id%2)+energy;
  pInfo.tally.cnt[id] = 0;

  pInfo.nInfo.rndState[id] = state;
}

__global__ void transport(MemStruct Info){
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  int nid = Info.nInfo.id[id];

}
#if defined(__TRACK)
__global__ void history(int numIso, multipole isotope, CMPTYPE* devicearray, MemStruct Info, unsigned num_src, unsigned devstep){
#else
__global__ void history(int numIso, multipole isotope, MemStruct Info, unsigned num_src, unsigned devstep){
#endif
  //try others when real simulation structure becomes clear
  int idl = threadIdx.x;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  //printf("i'm thread %d\n",id);
  int nid = Info.nInfo.id[id];
  unsigned live;
  unsigned isotopeID=Info.nInfo.isotope[nid];
  extern __shared__ unsigned blockTerminated[];
  CMPTYPE localenergy;
  CMPTYPE rnd;
  CMPTYPE sigT, sigA, sigF;
  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.nInfo.rndState[nid];

  localenergy = Info.nInfo.energy[nid];
  live = 1u;
  //while(live){
  //for (istep = 0; istep < devstep; istep++){
    rnd = curand_uniform(&localState);
#if defined(__SAMPLE)
    isotope.xs_eval_fast(localenergy + 
		      curand_normal(&localState)*sqrt(300.0*KB)*sqrt(0.5)/mp_para.dev_doubles[SQRTAWR], 
		      sigT, sigA, sigF);
#else
    isotope.xs_eval_fast(isotopeID,localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
    //isotope.xs_eval_fast(1-isotopeID, localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
#endif
#if defined(__TRACK)
    unsigned lies = gridDim.x*blockDim.x;
    live = Info.tally.cnt[nid] + cnt;
    live = live*(live<lies) + lies*(live>=lies); 
    if(0==id){
      devicearray[4*live  ] = localenergy;
      devicearray[4*live+1] = sigF;  
    }
    if(2==id){
      devicearray[4*live+2] = localenergy;
      devicearray[4*live+3] = sigF;
    }
#endif

    localenergy = localenergy * rnd;
    live = (localenergy > 1.0);
    isotopeID = rnd<0.5; //id%2;//0;//an example law to change isotopeID 
    /*So far, energy is the only state*/
    localenergy = localenergy*live + STARTENE*(1u - live);
    //terminated += !live;
  //}
  //}
  blockTerminated[idl] = !live;
   __syncthreads();
  live = blockDim.x>>1;
  while(live){
    if(idl<live)
      blockTerminated[idl] += blockTerminated[idl+live];
    __syncthreads();
    live>>=1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    Info.block_terminated_neutrons[blockIdx.x] = blockTerminated[0];
  }
  
  /*Note: from now on, live does not indicate neutron but thread active */
  //blockActive[threadIdx.x] = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  //atomicAdd(Info.num_terminated_neutrons,!live);
  //Info.thread_active[id] =  blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  Info.nInfo.rndState[nid] = localState; 
  Info.nInfo.energy[nid] = localenergy;
  Info.nInfo.isoenergy[id] = localenergy+isotopeID*MAXENERGY;
  Info.nInfo.isotope[nid] = isotopeID;
  Info.tally.cnt[nid] += 1; 

}


__global__ void remaining(int numIso,multipole isotope, CMPTYPE *devicearray, MemStruct Info){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  int nid=Info.nInfo.id[id];
  unsigned live = true;
  unsigned isotopeID=Info.nInfo.isotope[nid];
  CMPTYPE localenergy;
  CMPTYPE rnd;
  CMPTYPE sigT, sigA, sigF;
 
  /* Copy state to local memory for efficiency */
  curandState localState = Info.nInfo.rndState[nid];
  
#if defined(__PROCESS)
  localenergy = 1.0+19999.0/65536.0*id+0.181317676432466;
#else
  localenergy = Info.nInfo.energy[nid];
#endif
  unsigned cnt = 0u;
  unsigned terminated = 0u;
  live = 1u;
  while(live){
    rnd = curand_uniform(&localState);
#if defined(__SAMPLE)
    isotope.xs_eval_fast(localenergy + 
		      curand_normal(&localState)*sqrt(300.0*KB)*sqrt(0.5)/mp_para.dev_doubles[SQRTAWR], 
		      sigT, sigA, sigF);
#else
    isotope.xs_eval_fast(isotopeID, localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
    //isotope.xs_eval_fast(1-isotopeID, localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
#endif
#if defined(__TRACK)
    unsigned lies = gridDim.x*blockDim.x;
    live = Info.tally.cnt[nid] + cnt;
    live = live*(live<lies) + lies*(live>=lies); 
    if(0==id){
      devicearray[4*live  ] = localenergy;
      devicearray[4*live+1] = sigF;  
    }
    if(2==id){
      devicearray[4*live+2] = localenergy;
      devicearray[4*live+3] = sigF;
    }
#endif

#if !defined(__PROCESS)
    localenergy = localenergy * rnd;
    live = (localenergy > 1.0);
    isotopeID = rnd<0.5; //0;// an example law to change isotopeID
    cnt = cnt + 1;
    terminated += !live;
#else
    live = false;
#endif
  }
  /* Copy state back to global memory */
  atomicAdd(Info.num_terminated_neutrons,terminated);
  //Info.nInfo.rndState[nid] = localState;
  Info.tally.cnt[nid] += cnt;

#if !defined(__TRACK)
#if defined(__PROCESS)  
  devicearray[4 * nid] = localenergy ;
#else
  devicearray[4 * nid] = localenergy / rnd;
#endif
  devicearray[4 * nid + 1] = sigT;
  devicearray[4 * nid + 2] = sigA;
  devicearray[4 * nid + 3] = sigF;
#endif
}

__global__ void statistics(unsigned *threadcnt, unsigned* cnt){
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
  shared[idl] = threadcnt[id];
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
    cnt[blockIdx.x] += shared[0];
  }
  
}

/*
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
*/
