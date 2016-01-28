#include "manmemory.h"
__constant__ float wdspp[9];
//Simulation memory allocate and deallocate
void initialize_device(){
  cudaSetDevice(0);
  gpuErrchk(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceLmemResizeToMax));
}

void copymeans(int *h_cnt, int *batcnt, unsigned meshes, unsigned offset){
  for(int im=0;im<meshes;im++)
    batcnt[offset+im] = h_cnt[im];

}


void initialize_memory_data(MemStruct *DeviceMem, MemStruct *HostMem){
  (*HostMem).wdspp = (float*)malloc(sizeof(float)*9);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 9*sizeof(float)));
}
void copydata(MemStruct DeviceMem, MemStruct HostMem){
  gpuErrchk(cudaMemcpy(DeviceMem.wdspp,  HostMem.wdspp,   sizeof(float)*9, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(wdspp, DeviceMem.wdspp, 9*sizeof(float), 0, cudaMemcpyDeviceToDevice));
}
void release_memory_data(MemStruct DeviceMem, MemStruct HostMem){
  free(HostMem.wdspp);
  gpuErrchk(cudaFree(DeviceMem.wdspp));
}
//==============================================================================
//=====memory_converge() allocates memory for the phase of source convergence===
//1. memory_converge() need not have tally arrays, which exist current for debug
//==============================================================================
void allocate_memory_converge(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned num_seg){
  unsigned gridsize,banksize;
  gridsize = gridx*blockx;
  banksize = gridx*blockx*num_seg;

#if defined(__TALLY)
#if defined(__FTALLY)||(__FTALLY2) //Fission source tally
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  banksize*3*sizeof(int)));
  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*banksize);
  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  memset((*HostMem).batcnt, 0, sizeof(CMPTYPE)*numbins);
#endif
#if defined(__CTALLY) //Collision density (flux) tally
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(CMPTYPE)));
  (*HostMem).spectrum = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);  

  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt), numbins*sizeof(CMPTYPE)));
#if defined(__CTALLY2)
  (*HostMem).batcnt2    = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt2),numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).cnt2_t),numbins*gridsize*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt2), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt2,0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif

#endif

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       banksize*sizeof(int)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, banksize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));

#if defined(__1D)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),banksize*sizeof(float)));
#endif 

#if defined(__3D)
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*banksize);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),3*banksize*sizeof(float)));
#endif

  return;
}


void release_memory_converge(MemStruct DeviceMem, MemStruct HostMem){
#if defined(__TALLY)
#if defined(__FTALLY)||(__FTALLY2) //Fission source tally
  gpuErrchk(cudaFree(DeviceMem.nInfo.imat));
  free(HostMem.nInfo.live);
  free(HostMem.batcnt);
#endif
#if defined(__CTALLY) //Collision density (flux) tally
  gpuErrchk(cudaFree(DeviceMem.spectrum));
  free(HostMem.spectrum);
  free(HostMem.spectrum);  

  free(HostMem.batcnt);
  gpuErrchk(cudaFree(DeviceMem.batcnt));
#if defined(__CTALLY2)
  free(HostMem.batcnt2);
  gpuErrchk(cudaFree(DeviceMem.batcnt2));
  gpuErrchk(cudaFree(DeviceMem.cnt2_t));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt2));
#endif
  gpuErrchk(cudaFree(DeviceMem.block_spectrum));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt));
#endif
#endif
  gpuErrchk(cudaFree(DeviceMem.nInfo.live));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));

#if defined(__1D)
  gpuErrchk(cudaFree(DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaFree(DeviceMem).nInfo.pos_y),banksize*sizeof(float)));
#endif 

#if defined(__3D)
  free(HostMem.nInfo.pos_x);
  free(HostMem.nInfo.pos_y);
  free(HostMem.nInfo.pos_z);
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_z));
#endif
  return;
}


//==============================================================================
//====================initialize memory for  delayed bank ======================
//==============================================================================
void initialize_memory_bank(MemStruct *HostMem, unsigned banksize){
  (*HostMem).bank.x = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).bank.y = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).bank.z = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).bank.available = (unsigned*)malloc(sizeof(unsigned)*banksize);

  (*HostMem).bank.size   = (unsigned*)malloc(sizeof(unsigned));
  (*HostMem).bank.cursor_end       = (unsigned*)malloc(sizeof(unsigned));

  memset((*HostMem).bank.cursor_available, 0, sizeof(unsigned)*banksize);
  ((*HostMem).bank.size)[0]    = banksize;
  ((*HostMem).bank.cursor_end)[0] = 0; 
}

void release_memory_bank(MemStruct HostMem){
  free(HostMem.bank.x);
  free(HostMem.bank.y);
  free(HostMem.bank.z);

  free(HostMem.bank.size);
  free(HostMem.bank.cursor_end);
  free(HostMem.bank.cursor_available);
}


//==============================================================================
//=====memory_active() allocates memory for the phase of source convergence===
//1. memory_active() takes num_seg as parameter, which varies with phases ====
//2. memory_active() have tally array
//==============================================================================
void allocate_memory_active(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned num_seg){
  unsigned gridsize,banksize;
  gridsize = gridx*blockx;
  banksize = gridx*blockx*num_seg;

#if defined(__TALLY)
#if defined(__FTALLY)||(__FTALLY2) //Fission source tally
#if defined(__FTALLY)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  banksize*3*sizeof(int)));
  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*banksize);
#else //(__FTALLY2)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  gridsize*3*sizeof(int)));
  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*gridsize);
#endif
  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  memset((*HostMem).batcnt, 0, sizeof(CMPTYPE)*numbins);
#endif
#if defined(__CTALLY) //Collision density (flux) tally
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(CMPTYPE)));
  (*HostMem).spectrum = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);  

  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt), numbins*sizeof(CMPTYPE)));
#if defined(__CTALLY2)
  (*HostMem).batcnt2    = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt2),numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).cnt2_t),numbins*gridsize*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt2), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt2,0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
#endif

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       banksize*sizeof(int)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, banksize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));

#if defined(__1D)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),banksize*sizeof(float)));
#endif 

#if defined(__3D)
#if defined(__FTALLY2)
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*gridsize);
#else
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*banksize);
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),3*banksize*sizeof(float)));
#endif

  return;
}


void release_memory_active(MemStruct DeviceMem, MemStruct HostMem){
#if defined(__TALLY)
#if defined(__FTALLY)||(__FTALLY2) //Fission source tally
  gpuErrchk(cudaFree(DeviceMem.nInfo.imat));
  free(HostMem.nInfo.live);
  free(HostMem.batcnt);
#endif
#if defined(__CTALLY) //Collision density (flux) tally
  gpuErrchk(cudaFree(DeviceMem.spectrum));
  free(HostMem.spectrum);  

  free(HostMem.batcnt);
  gpuErrchk(cudaFree(DeviceMem.batcnt));
#if defined(__CTALLY2)
  free(HostMem.batcnt2);
  gpuErrchk(cudaFree(DeviceMem.batcnt2));
  gpuErrchk(cudaFree(DeviceMem.cnt2_t));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt2));
#endif
  gpuErrchk(cudaFree(DeviceMem.block_spectrum));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt));
#endif
#endif

  gpuErrchk(cudaFree(DeviceMem.nInfo.live));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));

#if defined(__1D)
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
#endif 

#if defined(__3D)
  free(HostMem.nInfo.pos_x);
  free(HostMem.nInfo.pos_y);
  free(HostMem.nInfo.pos_z);
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_z));
#endif

  return;
}

void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned nbat,unsigned ubat){
  unsigned gridsize,banksize;
  gridsize = gridx*blockx;
  banksize = gridx*blockx*ubat;
  //for __TALLY, ubat is used as tranfac

#if defined(__TALLY)
#if defined(__FTALLY)||(__FTALLY2) //Fission source tally
#if defined(__FTALLY)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  banksize*3*sizeof(int)));
  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*banksize);
#else //(__FTALLY2)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  gridsize*3*sizeof(int)));
  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*gridsize);
  printf("live allocated:%d\n",gridsize);
#endif
  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  memset((*HostMem).batcnt, 0, sizeof(CMPTYPE)*numbins);
#endif
#if defined(__CTALLY) //Collision density (flux) tally
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(CMPTYPE)));
  (*HostMem).spectrum = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);  

  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt), numbins*sizeof(CMPTYPE)));
#if defined(__CTALLY2)
  (*HostMem).batcnt2    = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt2),numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).cnt2_t),numbins*gridsize*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt2), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt2,0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
#if defined(__PROCESS)
  (*HostMem).batchmeans = (double*)malloc(sizeof(double)*nbat*numbins);
  (*HostMem).accmeans   = (double*)malloc(sizeof(double)*(nbat-ubat)*numbins);
#endif
#endif
#if defined(__TRAN)&&defined(__TALLY)
  memset((*HostMem).batcnt, 0, sizeof(CMPTYPE)*numbins);
#endif
  (*HostMem).wdspp = (float*)malloc(sizeof(float)*9);



  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 9*sizeof(float)));


  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       banksize*sizeof(int)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, banksize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));

#if defined(__WASTE)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       banksize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigT),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigA),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigF),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),banksize*sizeof(CMPTYPE)));
#endif

#if defined(__1D)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),banksize*sizeof(float)));
#endif 

#if defined(__3D)
#if defined(__FTALLY2)
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*gridsize);
#else
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*banksize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*banksize);
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),3*banksize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),3*banksize*sizeof(float)));
#endif

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_terminated_neutrons), sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).num_terminated_neutrons, 0, sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_live_neutrons), sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridx));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;


  return;
}

void resettally(CMPTYPE *cnt, unsigned totbins){
  gpuErrchk(cudaMemset(cnt, 0, totbins*sizeof(CMPTYPE)));
}

void resettally(int *cnt, unsigned totbins){
  gpuErrchk(cudaMemset(cnt, 0, totbins*sizeof(int)));
}

void release_memory(MemStruct DeviceMem, MemStruct HostMem){
#if defined(__TALLY)
#if defined(__MTALLY)||(__FTALLY)||(__FTALLY2)
  gpuErrchk(cudaFree(DeviceMem.nInfo.imat));
#else
  free(HostMem.spectrum);
  gpuErrchk(cudaFree(DeviceMem.batcnt));
#if defined(__CTALLY2)
  gpuErrchk(cudaFree(DeviceMem.batcnt2));
  free(HostMem.batcnt2);
  gpuErrchk(cudaFree(DeviceMem.cnt2_t));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt2));
#endif
  gpuErrchk(cudaFree(DeviceMem.spectrum));
  gpuErrchk(cudaFree(DeviceMem.block_spectrum));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt));
#endif
  free(HostMem.batcnt);
#if defined(__PROCESS)
  free(HostMem.batchmeans);
  free(HostMem.accmeans);
#endif
#endif
  free(HostMem.nInfo.live);

  gpuErrchk(cudaFree(DeviceMem.nInfo.live));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));

#if defined(__WASTE)
  gpuErrchk(cudaFree(DeviceMem.nInfo.id));
  gpuErrchk(cudaFree(DeviceMem.nInfo.energy));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigT));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigA));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigF));
  gpuErrchk(cudaFree(DeviceMem.nInfo.isoenergy));
#endif

#if defined(__1D)
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
#endif

#if defined(__3D)
  free(HostMem.nInfo.pos_x);
  free(HostMem.nInfo.pos_y);
  free(HostMem.nInfo.pos_z);
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_z));
#endif

  gpuErrchk(cudaFree(DeviceMem.num_terminated_neutrons));
  gpuErrchk(cudaFree(DeviceMem.num_live_neutrons));
  gpuErrchk(cudaFree(DeviceMem.block_terminated_neutrons));
  gpuErrchk(cudaFreeHost(HostMem.num_terminated_neutrons));


  return;
}

