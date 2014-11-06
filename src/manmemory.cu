#include "manmemory.h"
__constant__ unsigned spectrumbins[NUM_BINS+1];

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

void assign_tallybins(double *h_tallybins, double **d_tallybins){
  gpuErrchk(cudaMalloc((void**)(d_tallybins), (NUM_BINS+1)*sizeof(double)));
  gpuErrchk(cudaMemcpy(*d_tallybins,h_tallybins,(NUM_BINS+1)*sizeof(double),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(spectrumbins, *d_tallybins, (NUM_BINS+1)*sizeof(double), 0, cudaMemcpyDeviceToDevice));
}

void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned **h_blockcnt, unsigned** d_blockcnt, unsigned gridx, unsigned blockx ){
  unsigned gridsize;
  gridsize = gridx*blockx;

  *h_blockcnt      = (unsigned*)malloc(gridx*sizeof(unsigned));

  gpuErrchk(cudaMalloc((void**)(d_blockcnt), gridx*sizeof(unsigned int)));
  gpuErrchk(cudaMemset(*d_blockcnt, 0, gridx*sizeof(unsigned int)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.id),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.live),       gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.rndState), gridsize*sizeof(curandState)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.energy),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigT),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigA),   gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.sigF),   gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.isoenergy),gridsize*sizeof(CMPTYPE)));

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.imat),  gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).nInfo.pos_x),gridsize*sizeof(float)));
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

  gpuErrchk(cudaMalloc((void**)&((*DeviceMem).tally.cnt), gridsize*sizeof(unsigned)));
  gpuErrchk(cudaMemset((*DeviceMem).tally.cnt, 0, gridsize*sizeof(unsigned)));  

  return;
}

void release_memory(MemStruct DeviceMem, MemStruct HostMem, unsigned *h_blockcnt, unsigned* d_blockcnt, double* d_tallybins ){
  free(h_blockcnt);
  gpuErrchk(cudaFree(d_blockcnt));
  gpuErrchk(cudaFree(d_tallybins));

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

