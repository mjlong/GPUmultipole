#include "manmemory.h"
__constant__ float spectrumbins[NUM_BINS+1];

#if defined (__FOURIERW)
#if defined(__ALLCPU)
#include "fourierw.hh"
#else
#include "fourierw.h"
__constant__ CMPTYPE a[M+1];
__constant__ CMPTYPE b[M+1];
#endif
#endif

#if defined (__QUICKW)
#if defined (__ALLCPU)
#include "QuickW.hh"
#else
#include "QuickW.h"
#if defined (__QUICKWC) || defined(__INTERPEXP)
__constant__ CMPTYPE2 constwtable[LENGTH*LENGTH];
#endif
#endif
#endif



//Simulation memory allocate and deallocate
void initialize_device(){
  cudaSetDevice(0);
  gpuErrchk(cudaSetDeviceFlags(cudaDeviceMapHost | cudaDeviceLmemResizeToMax));
}

#if defined(__FOURIERW)
#if defined(__ALLCPU)
void fill_wtables(CMPTYPE** da, CMPTYPE** db){
  *da = (CMPTYPE*)malloc(sizeof(CMPTYPE)*(M+1));
  *db = (CMPTYPE*)malloc(sizeof(CMPTYPE)*(M+1)); 
  fill_a(*da,*db); 
}
void release_wtables(CMPTYPE* da, CMPTYPE* db){
  free(da);
  free(db);
}
#else
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
#if defined(__ALLCPU)
void fill_wtables(CPUComplex<CMPTYPE>** wtable){
  *wtable = (CPUComplex<CMPTYPE>*)malloc(LENGTH*LENGTH * 2 * sizeof(CMPTYPE));
  fill_w_tabulated(*wtable);
}
void release_wtables(CPUComplex<CMPTYPE>* wtable){
  free(wtable);
}
#else
void fill_wtables(CComplex<CMPTYPE>** wtable){
  gpuErrchk(cudaMalloc((void**)wtable, LENGTH*LENGTH * 2 * sizeof(CMPTYPE)));
  fill_w_tabulated<<<LENGTH,LENGTH>>>(*wtable);
}
void release_wtables(CComplex<CMPTYPE>* wtable){
  gpuErrchk(cudaFree(wtable));
}
#endif
#endif

void initialize_memory(MemStruct *HostMem, unsigned numbins){
  (*HostMem).spectrum = (unsigned*)malloc(sizeof(unsigned)*numbins);  
  memset((*HostMem).spectrum, 0, sizeof(unsigned)*numbins);
  return;
}

void release_memory(MemStruct HostMem){
  free(HostMem.spectrum);
  return;
}

#if defined(__XS_GPU)
void allocate_buffer(unsigned maxiso, unsigned** iS_d, 
                     CMPTYPE** sigTs_h, CMPTYPE** sigAs_h, CMPTYPE** sigFs_h,
                     CMPTYPE** sigTs_d, CMPTYPE** sigAs_d, CMPTYPE** sigFs_d){
  *sigTs_h = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
  *sigAs_h = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
  *sigFs_h = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
  gpuErrchk(cudaMalloc((void**)iS_d, maxiso*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)sigTs_d, maxiso*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)sigAs_d, maxiso*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)sigFs_d, maxiso*sizeof(CMPTYPE)));
}

void release_buffer(unsigned* iS_d, 
                    CMPTYPE* sigTs_h, CMPTYPE* sigAs_h, CMPTYPE* sigFs_h, 
                    CMPTYPE* sigTs_d, CMPTYPE* sigAs_d, CMPTYPE* sigFs_d){
  free(sigTs_h);
  free(sigAs_h);
  free(sigFs_h);

  gpuErrchk(cudaFree(iS_d));
  gpuErrchk(cudaFree(sigTs_d));
  gpuErrchk(cudaFree(sigAs_d));
  gpuErrchk(cudaFree(sigFs_d));
}
#endif

#if defined(__W__GPU)
#include "multipole_data.h"
void allocate_zwarray(CPUComplex<CMPTYPE>** z_h, CComplex<CMPTYPE>** z_d, 
                      CPUComplex<CMPTYPE>** w_h, CComplex<CMPTYPE>** w_d, 
                      unsigned numiso, struct multipoledata *iso){
//find maxwindow
  int maxwindow = 1;
  for(int i=0;i<numiso;i++){
    for(int j=0;j<iso[i].windows;j++){
      int window = iso[i].w_end[j]-iso[i].w_start[j]+1;
      if( window > maxwindow)
        maxwindow = window;
    }
    printf("iso[%d].maxwindow=%d\n",i,maxwindow);
  }
  printf("maxwindow=%d\n",maxwindow);
  *z_h = (CPUComplex<CMPTYPE>*)malloc(sizeof(CMPTYPE)*2*maxwindow);
  *w_h = (CPUComplex<CMPTYPE>*)malloc(sizeof(CMPTYPE)*2*maxwindow);
  gpuErrchk(cudaMalloc((void**)z_d, maxwindow*sizeof(CMPTYPE)*2));
  gpuErrchk(cudaMalloc((void**)w_d, maxwindow*sizeof(CMPTYPE)*2));

}
void release_zwarray(CPUComplex<CMPTYPE>* z_h, CComplex<CMPTYPE>* z_d, 
                     CPUComplex<CMPTYPE>* w_h, CComplex<CMPTYPE>* w_d){
  gpuErrchk(cudaFree(z_d));
  gpuErrchk(cudaFree(w_d));
  free(z_h);
  free(w_h);
}
#endif
