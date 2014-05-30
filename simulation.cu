#include "simulation.h"

__global__ void initialize(neutronInfo Info, double energy){
  //  int id = ((blockDim.x*blockDim.y*blockDim.z)*(blockIdx.y*gridDim.x+blockIdx.x)+(blockDim.x*blockDim.y)*threadIdx.z+blockDim.x*threadIdx.y+threadIdx.x);//THREADID;
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  Info.energy[id] = energy; //id+1.0; //(id + 1)*1.63*energy*0.001;// 

}

__global__ void history(multipole U238, double *devicearray, struct neutronInfo Info){
  //__global__ void history(multipole U238, struct neutronInfo Info){
  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = blockDim.x * blockIdx.x + threadIdx.x;//THREADID;

  bool live=true;
  double localenergy;
  double rnd;
  double sigT, sigA, sigF;
  struct pointers sharedptr;
  extern __shared__ float shared[];
  //size of shared[] is given as 3rd parameter while launching the kernel
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &Info.rndState[id]);

  /* Copy state to local memory for efficiency */ 
  curandState localState = Info.rndState[id];

  localenergy = Info.energy[id];
  unsigned cnt = 0;
  unsigned idl = threadIdx.x;
  unsigned idb = blockIdx.x;
  unsigned blocksize = blockDim.x;
  /*
    shift shared memory for double twophi[MAXNUML] and complex sigT_factor[MAXNUML]
  */
  //TODO: tailor to accomodate more than two isotopes
  sharedptr.blockbase   = Info.share.blockbase;
  sharedptr.sigT_factor = (CComplex*)(shared + blocksize) + idl;
  sharedptr.w_start     = (unsigned*)(shared + blocksize + (blocksize<<2)*Info.share.numL);
  sharedptr.w_end       = sharedptr.w_start + Info.share.windows;

  cnt = idl;
  /*
  while(cnt<Info.share.windows){
    sharedptr.w_start[cnt] = U238.w_start[cnt];
    sharedptr.w_end[cnt]   = U238.w_end[cnt];
    cnt += blocksize;
  }*/
  //__syncthreads();
  cnt = 0;

  sharedptr.dev_doubles =   (double*)(sharedptr.w_end + Info.share.windows);
  if(idl<3)
    sharedptr.dev_doubles[idl]  = U238.dev_doubles[idl];
  //__syncthreads();

  sharedptr.pseudo_rho  =   sharedptr.dev_doubles + 3;
  if(idl<Info.share.numL)
    sharedptr.pseudo_rho[idl]   = U238.pseudo_rho[idl];
  //__syncthreads();

  sharedptr.dev_integers =  (unsigned*)(sharedptr.pseudo_rho + Info.share.numL);
  if(idl<4)
    sharedptr.dev_integers[idl] = U238.dev_integers[idl];
  __syncthreads();
  
  while(live){
    rnd = curand_uniform(&localState);
    U238.xs_eval_fast(localenergy, sqrt(300.0*KB), sigT, sigA, sigF, sharedptr);
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
  unsigned *tally = (unsigned*)(shared);
  int i;
  tally[idl] = cnt;
  __syncthreads();
  i = blocksize>>1;
  while(i){
    if(idl<i)
      tally[idl] += tally[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl){
    //reduction scheme depends on tally type
    //following is to count moderation times
    Info.ntally.cnt[idb] = tally[0];
  }

}


