#include "manmemory.h"
__constant__ float wdspp[9];
//Simulation memory allocate and deallocate
void initialize_device(){
  cudaSetDevice(0);
  gpuErrchk(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceLmemResizeToMax));
}

void copymeans(int *h_cnt, int *batcnt, unsigned meshes, unsigned offset){
  for(int im=0;im<meshes;im++)
#if defined(__TRAN)
    batcnt[offset+im] += h_cnt[im];
#else
    batcnt[offset+im] = h_cnt[im];
#endif

}

void copydata(MemStruct DeviceMem, MemStruct HostMem){
  gpuErrchk(cudaMemcpy(DeviceMem.wdspp,  HostMem.wdspp,   sizeof(float)*9, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(wdspp, DeviceMem.wdspp, 9*sizeof(float), 0, cudaMemcpyDeviceToDevice));
}
void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned nbat,unsigned ubat){
  unsigned gridsize;
#if defined(__TRAN)
  (*HostMem).initial_delayed = (int*)malloc(sizeof(int)*nbat);  
  gridsize = gridx*blockx*ubat;
  //for __TALLY, ubat is used as tranfac

  (*HostMem).nInfo.d_pos_x = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.d_pos_y = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.d_pos_z = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.d_time  = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.d_igen  = (int*)malloc(sizeof(int)*gridsize);
  (*HostMem).nInfo.d_nu    = (int*)malloc(sizeof(int)*gridsize);
#else
  gridsize = gridx*blockx;
#endif

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


  (*HostMem).nInfo.live  = (int*)malloc(sizeof(int)*gridsize);

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 9*sizeof(float)));


  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.live, 0, gridsize*sizeof(int)));  

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));

#if defined(__WASTE)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigT),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigA),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigF),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  gridsize*sizeof(int)));
#endif

#if defined(__SCATTERPLOT)
  (*HostMem).nInfo.energy  = (CMPTYPE*)malloc(sizeof(CMPTYPE)*gridsize);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   gridsize*sizeof(CMPTYPE))); //use as initial z position for plot
#endif

#if defined(__1D)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),gridsize*sizeof(float)));
#endif 

#if defined(__3D)
  (*HostMem).nInfo.pos_x = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_y = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.pos_z = (float*)malloc(sizeof(float)*gridsize);
#if defined(__TRAN)
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),gridsize*sizeof(float)));

  (*HostMem).nInfo.dir_polar = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.dir_azimu = (float*)malloc(sizeof(float)*gridsize);
  (*HostMem).nInfo.d_closest = (float*)malloc(sizeof(float)*gridsize);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_polar),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_azimu),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.d_closest ),gridsize*sizeof(float)));
  gpuErrchk(cudaMemset((*DeviceMem).nInfo.d_closest, 0, gridsize*sizeof(float)));  //use as time
#else
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),3*gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),3*gridsize*sizeof(float)));
#endif
#endif

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_terminated_neutrons), sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).num_terminated_neutrons, 0, sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_live_neutrons), sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridx));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;


  return;
}

void resettally(int *cnt, unsigned totbins){
  gpuErrchk(cudaMemset(cnt, 0, totbins*sizeof(int)));}

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

#if defined(__TRAN)
  free(HostMem.nInfo.d_pos_x);
  free(HostMem.nInfo.d_pos_y);
  free(HostMem.nInfo.d_pos_z);
  free(HostMem.nInfo.d_time );
  free(HostMem.nInfo.d_nu   );
  free(HostMem.nInfo.d_igen );

  free(HostMem.initial_delayed);
  free(HostMem.nInfo.dir_polar);
  free(HostMem.nInfo.dir_azimu);
  free(HostMem.nInfo.d_closest);
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_polar));
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_azimu));
  gpuErrchk(cudaFree(DeviceMem.nInfo.d_closest));
#endif
#endif

  gpuErrchk(cudaFree(DeviceMem.num_terminated_neutrons));
  gpuErrchk(cudaFree(DeviceMem.num_live_neutrons));
  gpuErrchk(cudaFree(DeviceMem.block_terminated_neutrons));
  gpuErrchk(cudaFreeHost(HostMem.num_terminated_neutrons));


  return;
}

