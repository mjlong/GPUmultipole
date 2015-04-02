#include "simulation.h" 

__global__ void initialize(MemStruct pInfo,float width){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &(pInfo.nInfo.rndState[id]));

  neutron_sample(pInfo.nInfo, id,width);
  pInfo.nInfo.id[id] = id;
  pInfo.tally.cnt[id] = 0;
}

__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id,float width){
  nInfo.live[id] = 1u;
  curandState state = nInfo.rndState[id];
  //TODO: source sampling should take settings dependent on geometry
  nInfo.pos_x[id] = width*curand_uniform(&state);
  nInfo.pos_y[id] = 0.5f+0.00*curand_uniform(&state);
  nInfo.pos_z[id] = 0.5f+0.00*curand_uniform(&state);
  nInfo.dir_polar[id] = curand_uniform(&state)*2-1;
  nInfo.dir_azimu[id] = curand_uniform(&state)*PI*2;
  nInfo.energy[id] = STARTENE;
  nInfo.rndState[id] = state;
}


__global__ void history(MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned devstep,float width, float dx){
  //try others when real simulation structure becomes clear
  int idl = threadIdx.x;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  int nid = id;
  unsigned live;
  extern __shared__ unsigned blockTerminated[];

  CMPTYPE rnd;
  float x = DeviceMem.nInfo.pos_x[nid];
  int dir = 1-2*int(DeviceMem.nInfo.dir_polar[nid]<=0.5);
  /* Copy state to local memory for efficiency */ 
  curandState localState = DeviceMem.nInfo.rndState[nid];

  unsigned istep;
  //printf("[%2d],x=%.5f\n",id,DeviceMem.nInfo.pos_x[nid]);
  for(istep=0;istep<devstep;istep++){
    DeviceMem.tally.cnt[int(x/dx)*gridDim.x*blockDim.x+nid]+=1;
    rnd = curand_uniform(&localState);
    

    rnd = curand_uniform(&localState);

    DeviceMem.nInfo.live[nid] = live;  
    //energy can be updated efficiently here, live state is upated after sorting
    live = rnd<0.5;
    //terminated += !live;

  }
  blockTerminated[idl] = !live;
  
  /*Note: from now on, live does not indicate neutron but thread active */
  //blockActive[threadIdx.x] = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  //atomicAdd(Info.num_terminated_neutrons,!live);
  //Info.thread_active[id] =  blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  DeviceMem.nInfo.rndState[nid] = localState; 

  /*
  else{
    blockTerminated[idl] = active;//0;
    //those old unlive neutrons must not be counted again
    //so, 0 instead of !live is used 
    //it was incorrect, above senario forgot to count leak neutron as terminated
  }
  */
  //TODO: no need of such within block reduction for remaining()
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
    DeviceMem.block_terminated_neutrons[blockIdx.x] = blockTerminated[0];
  }

}

__global__ void reduce_sum_plus(unsigned *threadcnt, unsigned* cnt){
// reduce threadcnt[] to cnt[], cnt is updated by self increase
// this is used to count terminated neurtons
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

__global__ void reduce_sum_equal(unsigned* thread_active, unsigned* active){
// reduce thread_active to active, active is updated without history
// this is used to count number of "live" threads
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned idl = threadIdx.x;
  extern __shared__ unsigned shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  int i;
  shared[idl] = thread_active[id]; 
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
