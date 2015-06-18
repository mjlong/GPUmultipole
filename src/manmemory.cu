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

void copydata(MemStruct DeviceMem, MemStruct HostMem){
  gpuErrchk(cudaMemcpy(DeviceMem.wdspp,  HostMem.wdspp,   sizeof(float)*9, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(wdspp, DeviceMem.wdspp, 9*sizeof(float), 0, cudaMemcpyDeviceToDevice));
}

void delayed_memory(int nbat, int num_srcp, int csize,MemStruct* HostMem){
  memset((*HostMem).newly_delayed,   0, sizeof(int)*nbat);
  (*HostMem).nInfo.d_pos_x = (float*)malloc(sizeof(float)*csize);
  (*HostMem).nInfo.d_pos_y = (float*)malloc(sizeof(float)*csize);
  (*HostMem).nInfo.d_pos_z = (float*)malloc(sizeof(float)*csize);
  (*HostMem).nInfo.d_igen  = (int*)malloc(sizeof(int)*csize);
  (*HostMem).nInfo.d_nu    = (int*)malloc(sizeof(int)*csize);
  for(int i=0;i<csize;i++)
    ((*HostMem).nInfo.d_igen)[i] = -1;
  (*HostMem).nInfo.dbank_x = (float*)malloc(sizeof(float)*num_srcp);
  (*HostMem).nInfo.dbank_y = (float*)malloc(sizeof(float)*num_srcp);
  (*HostMem).nInfo.dbank_z = (float*)malloc(sizeof(float)*num_srcp);
}

void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned nbat,unsigned ubat){
  unsigned gridsize,banksize;
  gridsize = gridx*blockx;
  banksize = gridx*blockx*ubat;
  //for __TALLY, ubat is used as tranfac

  (*HostMem).newly_delayed = (int*)malloc(sizeof(int)*nbat);
#if defined(__TALLY)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(CMPTYPE)));
  (*HostMem).spectrum = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);  
  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt), numbins*sizeof(CMPTYPE)));
#if defined(__PROCESS)
  (*HostMem).batchmeans = (double*)malloc(sizeof(double)*nbat*numbins);
  (*HostMem).accmeans   = (double*)malloc(sizeof(double)*(nbat-ubat)*numbins);
#endif
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridx*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridx*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsize*sizeof(CMPTYPE)));  
#endif
#if defined(__TRAN)&&defined(__TALLY)
  memset((*HostMem).batcnt, 0, sizeof(CMPTYPE)*numbins);
#endif
  (*HostMem).wdspp = (float*)malloc(sizeof(float)*9);


  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*banksize);

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 9*sizeof(float)));


  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       banksize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       banksize*sizeof(int)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, banksize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), banksize*sizeof(curandState)));

#if defined(__WASTE)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigT),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigA),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigF),   banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),banksize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  banksize*sizeof(int)));
#endif

#if defined(__SCATTERPLOT)
  (*HostMem).nInfo.energy  = (CMPTYPE*)malloc(sizeof(CMPTYPE)*banksize);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   banksize*sizeof(CMPTYPE))); //use as initial z position for plot
#endif

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

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_terminated_neutrons), sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).num_terminated_neutrons, 0, sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_live_neutrons), sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridx));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;


  return;
}

void resettally(CMPTYPE *cnt, unsigned totbins){
  gpuErrchk(cudaMemset(cnt, 0, totbins*sizeof(CMPTYPE)));}

void release_memory(MemStruct DeviceMem, MemStruct HostMem){
  free(HostMem.nInfo.dbank_x);
  free(HostMem.nInfo.dbank_y);
  free(HostMem.nInfo.dbank_z);

  free(HostMem.nInfo.d_pos_x);
  free(HostMem.nInfo.d_pos_y);
  free(HostMem.nInfo.d_pos_z);
  free(HostMem.nInfo.d_nu   );
  free(HostMem.nInfo.d_igen );
  free(HostMem.newly_delayed);
#if defined(__TALLY)
  free(HostMem.spectrum);
#if defined(__PROCESS)
  free(HostMem.batchmeans);
  free(HostMem.accmeans);
#endif
  free(HostMem.batcnt);
  gpuErrchk(cudaFree(DeviceMem.batcnt));
  gpuErrchk(cudaFree(DeviceMem.spectrum));
  gpuErrchk(cudaFree(DeviceMem.block_spectrum));
  gpuErrchk(cudaFree(DeviceMem.tally.cnt));
#endif
  free(HostMem.nInfo.live);
  free(HostMem.wdspp);

  gpuErrchk(cudaFree(DeviceMem.wdspp));

  gpuErrchk(cudaFree(DeviceMem.nInfo.id));
  gpuErrchk(cudaFree(DeviceMem.nInfo.live));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));

#if defined(__WASTE)
  gpuErrchk(cudaFree(DeviceMem.nInfo.energy));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigT));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigA));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigF));
  gpuErrchk(cudaFree(DeviceMem.nInfo.isoenergy));
  gpuErrchk(cudaFree(DeviceMem.nInfo.imat));
#endif

#if defined(__1D)
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
#endif

#if defined(__SCATTERPLOT)
  gpuErrchk(cudaFree(DeviceMem.nInfo.energy));
  free(HostMem.nInfo.energy);
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

