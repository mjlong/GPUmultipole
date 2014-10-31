#include "material.h"

material::material(struct matdata *pmat, unsigned numIso){
  unsigned numMat = pmat->numMat;
  gpuErrchk(cudaMalloc((void**)&offsets,sizeof(unsigned)*(numMat+1))); 
  gpuErrchk(cudaMalloc((void**)&N_tot, sizeof(float)*numMat));
  gpuErrchk(cudaMalloc((void**)&densities, sizeof(float)*numIso));
  gpuErrchk(cudaMalloc((void**)&isotopes, sizeof(unsigned)*numIso));

  gpuErrchk(cudaMemcpy(offsets, pmat->offsets, sizeof(unsigned)*numMat, cudaMemcpyHostToDevice)); 
  shiftoffset<<<1,numMat>>>(offsets,numMat-1-(pmat->offsets[numMat-1]));
  gpuErrchk(cudaMemcpy(N_tot, pmat->N_tot, sizeof(float)*numMat, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(densities, pmat->densities, sizeof(float)*numIso, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(isotopes, pmat->isotopes, sizeof(unsigned)*numIso, cudaMemcpyHostToDevice));
}

material::~material(){
}

void material::release_pointer(){
  gpuErrchk(cudaFree(offsets));
  gpuErrchk(cudaFree(N_tot));
  gpuErrchk(cudaFree(densities));
  gpuErrchk(cudaFree(isotopes));
}

__global__ void shiftoffset(unsigned* offsets, unsigned shift)
{
  int id = threadIdx.x;
  offsets[id]+=shift;
}
