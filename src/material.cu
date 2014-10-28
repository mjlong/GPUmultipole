#include "material.h"

material::material(struct matdata *pmat, unsigned numIso){
  gpuErrchk(cudaMalloc((void**)&offsets,sizeof(unsigned)*pmat->numMat)); 
  gpuErrchk(cudaMalloc((void**)&N_tot, sizeof(float)*pmat->numMat));
  gpuErrchk(cudaMalloc((void**)&densities, sizeof(float)*numIso));
  gpuErrchk(cudaMalloc((void**)&isotopes, sizeof(unsigned)*numIso));

  gpuErrchk(cudaMemcpy(offsets, pmat->offsets, sizeof(unsigned)*pmat->numMat, cudaMemcpyHostToDevice)); 
  gpuErrchk(cudaMemcpy(N_tot, pmat->N_tot, sizeof(float)*pmat->numMat, cudaMemcpyHostToDevice));
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
