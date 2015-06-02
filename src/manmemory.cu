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
void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned ubat,unsigned gridr,unsigned blockr,unsigned ubatr){
  unsigned banksize,gridsizr;
  banksize = gridx*blockx*ubat;
  gridsizr = gridr*blockr;
  //for __TALLY, ubat is used as tranfac

#if defined(__TALLY)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(CMPTYPE)));
  (*HostMem).spectrum = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);  
  (*HostMem).batcnt     = (CMPTYPE*)malloc(sizeof(CMPTYPE)*numbins);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).batcnt), numbins*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridr*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridr*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsizr*numbins*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsizr*sizeof(CMPTYPE)));  
#endif
  (*HostMem).wdspp = (float*)malloc(sizeof(float)*9);


  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*banksize);

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 9*sizeof(float)));


  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       banksize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       banksize*sizeof(int)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, banksize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), banksize*sizeof(curandState)));


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

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridr));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;


  return;
}

void resettally(CMPTYPE *cnt, unsigned totbins){
  gpuErrchk(cudaMemset(cnt, 0, totbins*sizeof(CMPTYPE)));}

void release_memory(MemStruct DeviceMem, MemStruct HostMem){
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

