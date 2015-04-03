#include "manmemory.h"

//Simulation memory allocate and deallocate
void initialize_device(){
  cudaSetDevice(0);
  gpuErrchk(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceLmemResizeToMax));
}

void copymeans(unsigned *h_cnt, unsigned *acccnt, unsigned meshes, unsigned offset){
  for(int im=0;im<meshes;im++)
    acccnt[offset+im] = h_cnt[im];

}

void copydata(MemStruct DeviceMem, MemStruct HostMem){
  gpuErrchk(cudaMemcpy(DeviceMem.wdspp,  HostMem.wdspp,   sizeof(float)*5, cudaMemcpyHostToDevice));
}
void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned nbat){
  unsigned gridsize;
  gridsize = gridx*blockx;

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).spectrum), numbins*sizeof(unsigned int)));
  (*HostMem).spectrum = (unsigned*)malloc(sizeof(unsigned)*numbins);  
  (*HostMem).batchmeans = (float*)malloc(sizeof(float)*nbat*numbins);
  (*HostMem).accmeans   = (float*)malloc(sizeof(float)*nbat*numbins);
  (*HostMem).acccnt     = (unsigned*)malloc(sizeof(unsigned)*nbat*numbins);
  (*HostMem).wdspp = (float*)malloc(sizeof(float)*5);
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).wdspp), 5*sizeof(float)));
  
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_spectrum), numbins*gridx*sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).block_spectrum, 0, numbins*gridx*sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigT),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigA),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigF),   gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),3*gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_polar),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_azimu),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.d_closest ),gridsize*sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_terminated_neutrons), sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).num_terminated_neutrons, 0, sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_live_neutrons), sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridx));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*numbins*sizeof(unsigned)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, numbins*gridsize*sizeof(unsigned)));  

  return;
}

void release_memory(MemStruct DeviceMem, MemStruct HostMem){
  free(HostMem.spectrum);
  free(HostMem.batchmeans);
  free(HostMem.accmeans);
  free(HostMem.acccnt);
  free(HostMem.wdspp);

  gpuErrchk(cudaFree(DeviceMem.wdspp));
  gpuErrchk(cudaFree(DeviceMem.spectrum));
  gpuErrchk(cudaFree(DeviceMem.block_spectrum));

  gpuErrchk(cudaFree(DeviceMem.nInfo.id));
  gpuErrchk(cudaFree(DeviceMem.nInfo.live));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));
  gpuErrchk(cudaFree(DeviceMem.nInfo.energy));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigT));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigA));
  gpuErrchk(cudaFree(DeviceMem.nInfo.sigF));

  gpuErrchk(cudaFree(DeviceMem.nInfo.isoenergy));

  gpuErrchk(cudaFree(DeviceMem.nInfo.imat));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_z));
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_polar));
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_azimu));
  gpuErrchk(cudaFree(DeviceMem.nInfo.d_closest));

  gpuErrchk(cudaFree(DeviceMem.num_terminated_neutrons));
  gpuErrchk(cudaFree(DeviceMem.num_live_neutrons));
  gpuErrchk(cudaFree(DeviceMem.block_terminated_neutrons));
  gpuErrchk(cudaFreeHost(HostMem.num_terminated_neutrons));

  gpuErrchk(cudaFree(DeviceMem.tally.cnt));

  return;
}

