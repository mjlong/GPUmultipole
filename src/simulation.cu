#include "simulation.h" 
extern __constant__ float spectrumbins[];


__global__ void z2w_d(CComplex<CMPTYPE> *pz, CComplex<CMPTYPE> *pw,unsigned numz){
  unsigned id = blockDim.x*blockIdx.x+threadIdx.x;
  if(id<numz){
    pw[id] = Faddeeva::w(pz[id]);
  }
  //printf("w(%+.5e%+.5ei)=%+.5e%+.5ei\n",real(pz[id]),imag(pz[id]),real(pw[id]),imag(pw[id]));
}

__global__ void initialize(MemStruct pInfo){
  //int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &(pInfo.nInfo.rndState[id]));

  neutron_sample(pInfo.nInfo, id);
  pInfo.nInfo.id[id] = id;
  pInfo.tally.cnt[id] = 0;
}

__global__ void update_sort_key(MemStruct DeviceMem, material mat){
  unsigned id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned isoID = mat.isotopes[mat.offsets[DeviceMem.nInfo.imat[id]]+0];
                                            //matID
  DeviceMem.nInfo.isoenergy[id] = (MAXENERGY*isoID+DeviceMem.nInfo.energy[id])*DeviceMem.nInfo.live[id];
}

__global__ void transport(MemStruct DeviceMem, material mat,unsigned renew){
  int nid = (DeviceMem.nInfo.id[blockDim.x * blockIdx.x + threadIdx.x])%(gridDim.x*blockDim.x);
  unsigned live = DeviceMem.nInfo.live[nid];
  if(live){
    CMPTYPE sigT = DeviceMem.nInfo.sigT[nid];
    float s = -log(curand_uniform(&(DeviceMem.nInfo.rndState[nid])))/1E8;//(mat.N_tot[DeviceMem.nInfo.imat[nid]]*sigT);   
    float d = DeviceMem.nInfo.d_closest[nid];
    s = (d<s)*d+(d>=s)*s;
    float mu = DeviceMem.nInfo.dir_polar[nid];
    float phi= DeviceMem.nInfo.dir_azimu[nid];
    DeviceMem.nInfo.pos_x[nid]+=s*sqrt(1-mu*mu)*cos(phi);
    DeviceMem.nInfo.pos_y[nid]+=s*sqrt(1-mu*mu)*sin(phi);
    DeviceMem.nInfo.pos_z[nid]+=s*mu;
  }
  else if(renew){
    neutron_sample(DeviceMem.nInfo,nid);
    DeviceMem.nInfo.id[blockDim.x * blockIdx.x + threadIdx.x] += gridDim.x*blockDim.x;
  }
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


__global__ void resurrection(NeutronInfoStruct nInfo){
  //neutron energy has been set in an efficient way after each collison
  //only position and direction are sampled as neutron 
  unsigned nid = (nInfo.id[blockDim.x*blockIdx.x + threadIdx.x])%(gridDim.x*blockDim.x);
  unsigned live = nInfo.live[nid];
  if(!live)
    neutron_sample(nInfo,nid);
}
__device__ unsigned search_bin(CMPTYPE energy){
  for(int i=0;i<NUM_BINS;i++){
    if( (spectrumbins[i]>=energy)&&(spectrumbins[i+1]<energy) ) 
      return i;
  }
  return 0;
}
__global__ void history(material mat, multipole mp_para, MemStruct DeviceMem, unsigned num_src,unsigned active){
  //try others when real simulation structure becomes clear
  int idl = threadIdx.x;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  int nid = (DeviceMem.nInfo.id[id])%(gridDim.x*blockDim.x);
  unsigned live;
  extern __shared__ unsigned blockTerminated[];
  if(DeviceMem.nInfo.live[nid]){
  unsigned isotopeID;
  CMPTYPE localenergy;
  CMPTYPE rnd;
  CMPTYPE sigTsum, sigAsum, sigFsum, sigT, sigA, sigF;
  sigTsum=0;
  sigAsum=0;
  sigFsum=0;

  /* Copy state to local memory for efficiency */ 
  curandState localState = DeviceMem.nInfo.rndState[nid];

  localenergy = DeviceMem.nInfo.energy[nid];
  DeviceMem.tally.cnt[search_bin(localenergy)*gridDim.x*blockDim.x+nid]+=1;
  live = 1u;
  unsigned imat = DeviceMem.nInfo.imat[nid];
  for(isotopeID=mat.offsets[imat];isotopeID<mat.offsets[imat+1];isotopeID++ ){
    mp_para.xs_eval_fast(mat.isotopes[isotopeID],localenergy, sqrt(300.0*KB), sigT, sigA, sigF);
    sigTsum += sigT*mat.densities[isotopeID];
    sigAsum += sigA*mat.densities[isotopeID];
    sigFsum += sigF*mat.densities[isotopeID];
  }
  rnd = curand_uniform(&localState);

#if defined(__PRINTTRACK__)
  if(__PRINTTRACK__){
    printf("%7d[%5d],%3d,%+.7e, %+.7e, %+.7e, %.14e %.14e %.14e %.14e\n",
            DeviceMem.nInfo.id[id],id, DeviceMem.nInfo.imat[nid],
            DeviceMem.nInfo.pos_x[nid], DeviceMem.nInfo.pos_y[nid], DeviceMem.nInfo.pos_z[nid],
            localenergy, sigTsum,sigAsum,sigFsum); 
  }
#endif
  localenergy = localenergy * rnd;
  live = (localenergy > ENDENERG);
  DeviceMem.nInfo.live[nid] = live;  
  //energy can be updated efficiently here, live state is upated after sorting
  localenergy = localenergy*live + STARTENE*(1u - live);
  //terminated += !live;

  blockTerminated[idl] = !live;
  
  /*Note: from now on, live does not indicate neutron but thread active */
  //blockActive[threadIdx.x] = (((terminated*2)*blockDim.x*gridDim.x + atomicAdd(Info.num_terminated_neutrons, terminated)) < num_src);
  //atomicAdd(Info.num_terminated_neutrons,!live);
  //Info.thread_active[id] =  blockDim.x*gridDim.x + *Info.num_terminated_neutrons < num_src;
  /* Copy state back to global memory */ 
  DeviceMem.nInfo.rndState[nid] = localState; 
  DeviceMem.nInfo.energy[nid] = localenergy;
  DeviceMem.nInfo.sigT[nid]=sigTsum;
  DeviceMem.nInfo.sigA[nid]=sigAsum;
  DeviceMem.nInfo.sigF[nid]=sigFsum;
  }//end if live

  else{
    blockTerminated[idl] = active;//0;
    //those old unlive neutrons must not be counted again
    //so, 0 instead of !live is used 
    //it was incorrect, above senario forgot to count leak neutron as terminated
  }
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
