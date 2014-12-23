#include "manmemory.h"
__constant__ float spectrumbins[NUM_BINS+1];

#if defined (__FOURIERW)
#include "fourierw.h"
__constant__ CMPTYPE a[M+1];
__constant__ CMPTYPE b[M+1];
#endif

#if defined (__QUICKW)
#include "QuickW.h"
#endif

#if defined (__QUICKWC) || defined(__INTERPEXP)
__constant__ CMPTYPE2 constwtable[LENGTH*LENGTH];
#endif


//Simulation memory allocate and deallocate
void initialize_device(){
  cudaSetDevice(0);
  gpuErrchk(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceLmemResizeToMax));
}

#if defined(__FOURIERW)
void fill_wtables(CMPTYPE** da, CMPTYPE** db){
  gpuErrchk(cudaMalloc((void**)da, (M+1)*sizeof(CMPTYPE))); 
  gpuErrchk(cudaMalloc((void**)db, (M+1)*sizeof(CMPTYPE))); 
  fill_a<<<1,M+1>>>(*da,*db); 
  cudaMemcpyToSymbol(a, *da, M*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
  cudaMemcpyToSymbol(b, *db, M*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);

}
void release_wtables(CMPTYPE* da, CMPTYPE* db){
  gpuErrchk(cudaFree(da));
  gpuErrchk(cudaFree(db));
}
#endif

#if defined(__INTERPEXP)
void fill_wtables(CComplex<CMPTYPE>** exptable){
  gpuErrchk(cudaMalloc((void**)exptable, LENGTH*LENGTH * 2 * sizeof(CMPTYPE)));
  fill_exp_table<<<LENGTH,LENGTH>>>(*exptable);
  cudaMemcpyToSymbol(constwtable, *exptable, LENGTH*LENGTH*2*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
}
void release_wtables(CComplex<CMPTYPE>* exptable){
  gpuErrchk(cudaFree(exptable));
}
#endif

#if defined(__QUICKW)
void fill_wtables(CComplex<CMPTYPE>** wtable){
  gpuErrchk(cudaMalloc((void**)wtable, LENGTH*LENGTH * 2 * sizeof(CMPTYPE)));
  fill_w_tabulated<<<LENGTH,LENGTH>>>(*wtable);
#if defined(__QUICKWC)
  cudaMemcpyToSymbol(constwtable, *wtable, LENGTH*LENGTH*2*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
#endif
#if defined(__QUICKWT)
  bindwtable(*wtable);
#endif
}
void release_wtables(CComplex<CMPTYPE>* wtable){
  gpuErrchk(cudaFree(wtable));
}
#endif

void initialize_memory(MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx ){
  (*HostMem).spectrum = (unsigned*)malloc(sizeof(unsigned)*numbins);  
  (*HostMem).num_terminated_neutrons = (unsigned*)malloc(sizeof(unsigned));
  (*HostMem).num_terminated_neutrons[0] = 0u;
  return;
}

void release_memory(MemStruct HostMem){
  free(HostMem.spectrum);
  return;
}

