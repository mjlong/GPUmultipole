#include "manmemory.h"

//Simulation memory allocate and deallocate

void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, CMPTYPE** devicearray, CMPTYPE** hostarray, unsigned **cnt, unsigned** blockcnt, unsigned gridx, unsigned blockx ){
  unsigned gridsize;
  gridsize = gridx*blockx;

  *hostarray = (CMPTYPE*)malloc(4*gridsize*sizeof(CMPTYPE));
  *cnt      = (unsigned*)malloc(gridx*sizeof(unsigned));

  gpuErrchk(cudaMalloc((void**)devicearray, 4*gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset(*devicearray, 0, 4*gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)(blockcnt), gridx*sizeof(unsigned int)));
  gpuErrchk(cudaMemset(*blockcnt, 0, gridx*sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isotope),  gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.icell),  gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_y),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_z),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_polar),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.dir_azimu),gridsize*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.d_closest ),gridsize*sizeof(float)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).num_terminated_neutrons), sizeof(unsigned int)));
  gpuErrchk(cudaMemset((*DeviceMem).num_terminated_neutrons, 0, sizeof(unsigned)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).block_terminated_neutrons), sizeof(unsigned int)*gridx));
  gpuErrchk(cudaMallocHost((void**)&((*HostMem).num_terminated_neutrons), sizeof(unsigned int)));
  (*HostMem).num_terminated_neutrons[0] = 0u;

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, gridsize*sizeof(unsigned)));  


  return;
}

void release_memory(MemStruct DeviceMem, MemStruct HostMem, CMPTYPE* devicearray, CMPTYPE* hostarray, unsigned *cnt, unsigned* blockcnt ){
  free(hostarray);
  free(cnt);
  gpuErrchk(cudaFree(devicearray));
  gpuErrchk(cudaFree(blockcnt));

  gpuErrchk(cudaFree(DeviceMem.nInfo.id));
  gpuErrchk(cudaFree(DeviceMem.nInfo.rndState));
  gpuErrchk(cudaFree(DeviceMem.nInfo.energy));
  gpuErrchk(cudaFree(DeviceMem.nInfo.isotope));
  gpuErrchk(cudaFree(DeviceMem.nInfo.isoenergy));

  gpuErrchk(cudaFree(DeviceMem.nInfo.icell));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_x));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_y));
  gpuErrchk(cudaFree(DeviceMem.nInfo.pos_z));
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_polar));
  gpuErrchk(cudaFree(DeviceMem.nInfo.dir_azimu));
  gpuErrchk(cudaFree(DeviceMem.nInfo.d_closest));

  gpuErrchk(cudaFree(DeviceMem.num_terminated_neutrons));
  gpuErrchk(cudaFree(DeviceMem.block_terminated_neutrons));
  gpuErrchk(cudaFreeHost(HostMem.num_terminated_neutrons));

  gpuErrchk(cudaFree(DeviceMem.tally.cnt));
  return;
}

