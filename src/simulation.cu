#include "simulation.h" 
extern __constant__ float spectrumbins[];

__global__ void initialize(MemStruct pInfo){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &(pInfo.nInfo.rndState[id]));

  neutron_sample(pInfo.nInfo, id);
  pInfo.nInfo.id[id] = id;
  pInfo.tally.cnt[id] = 0;
}

__device__ void neutron_sample(NeutronInfoStruct nInfo, unsigned id){
  nInfo.live[id] = 1u;
  curandState state = nInfo.rndState[id];
//TODO: source sampling should take settings dependent on geometry
  nInfo.pos_x[id] = 0.5f+0.00*curand_uniform(&state);
  nInfo.pos_y[id] = 0.5f+0.00*curand_uniform(&state);
  nInfo.pos_z[id] = 0.5f+0.00*curand_uniform(&state);
  nInfo.dir_polar[id] = curand_uniform(&state)*2-1;
  nInfo.dir_azimu[id] = curand_uniform(&state)*PI*2;
  nInfo.energy[id] = STARTENE;
  nInfo.rndState[id] = state;
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
