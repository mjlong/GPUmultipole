#include "multipole.h"

/*multipole::multipole(){
  printf("Hello, i'm constructing\n");
}*/

#if defined(__QUICKWG)
multipole::multipole(struct multipoledata *data, int numIso, CComplex<CMPTYPE>* wtable){
#else 
multipole::multipole(struct multipoledata *data, int numIso){
#endif //Only __QUICWG needs a global wtable

  size_t size;
  int i;
  int * h_offset = (int*)malloc(sizeof(int)*numIso);
  int * h_size   = (int*)malloc(sizeof(int)*numIso);
  // allocate array of offsets
  size = sizeof(int)*numIso*NUMOFFS;
  gpuErrchk(cudaMalloc((void**)&offsets, size));
  /*
    allocate and assign integers
  */
  size = sizeof(unsigned);
  gpuErrchk(cudaMalloc((void**)&dev_integers, DEVINTS*size*numIso));
  for(i=0;i<numIso;i++){
    cudaMemcpy(dev_integers+i*DEVINTS+MODE,        &(data[i].mode), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_integers+i*DEVINTS+FITORDER,    &(data[i].fitorder), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_integers+i*DEVINTS+NUML,        &(data[i].numL), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_integers+i*DEVINTS+FISSIONABLE, &(data[i].fissionable), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_integers+i*DEVINTS+WINDOWS,     &(data[i].windows),size, cudaMemcpyHostToDevice);
  }
    cudaMalloc((void**)&dev_numIso,sizeof(int));
    cudaMemcpy(dev_numIso, &numIso, sizeof(int), cudaMemcpyHostToDevice);
  /*
    allocate and assign doubles
  */
  size = sizeof(CMPTYPE);
  cudaMalloc((void**)&dev_doubles,  DEVREALS*size*numIso);
  for(i=0;i<numIso;i++){
    cudaMemcpy(dev_doubles+i*DEVREALS+STARTE,  &(data[i].startE),  size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_doubles+i*DEVREALS+SPACING ,&(data[i].spacing), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_doubles+i*DEVREALS+SQRTAWR, &(data[i].sqrtAWR), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_doubles+i*DEVREALS+ENDE,    &(data[i].endE   ), size, cudaMemcpyHostToDevice);
  }

  /*
    allocate and assign arrays
  */
  // mpdata
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].length*(MP_RF+data[i].fissionable);
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=2*sizeof(CMPTYPE);
    size += h_size[i];
  }
    h_size[i] = data[i].length*(MP_RF+data[i].fissionable)*2*sizeof(CMPTYPE);
    size += h_size[i];
  cudaMalloc((void**)&mpdata, size);
  for(i=0;i<numIso;i++){
    cudaMemcpy(mpdata + h_offset[i], data[i].mpdata, h_size[i], cudaMemcpyHostToDevice);
  }
  cudaMemcpy(offsets+PMPDATA*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);

  // l_value
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].length;
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=sizeof(unsigned);
    size += h_size[i];
  }
    h_size[i] = data[i].length*sizeof(unsigned);
    size += h_size[i];
  cudaMalloc((void**)&l_value, size);
  for(i=0;i<numIso;i++){
    cudaMemcpy(l_value + h_offset[i], data[i].l_value, h_size[i], cudaMemcpyHostToDevice);
  }
  cudaMemcpy(offsets+PLVAL*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);

  // pseudo_rho
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].numL;
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=sizeof(CMPTYPE);
    size += h_size[i];
  }
    h_size[i] = data[i].numL*sizeof(CMPTYPE);
    size += h_size[i];
  cudaMalloc((void**)&pseudo_rho, size);
  for(i=0;i<numIso;i++){
    cudaMemcpy(pseudo_rho + h_offset[i], data[i].pseudo_rho, h_size[i], cudaMemcpyHostToDevice);
  }
  cudaMemcpy(offsets+PPRHO*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);

  // w_start and w_end
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].windows;
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=sizeof(int);
    size += h_size[i];
  }
    h_size[i] = data[i].windows*sizeof(int);
    size += h_size[i];
  cudaMalloc((void**)&w_start, size);
  cudaMalloc((void**)&w_end,   size);
  for(i=0;i<numIso;i++){
    cudaMemcpy(w_start + h_offset[i], data[i].w_start, h_size[i], cudaMemcpyHostToDevice);
    cudaMemcpy(w_end   + h_offset[i], data[i].w_end,   h_size[i], cudaMemcpyHostToDevice);
  }
  cudaMemcpy(offsets+PWIND*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);

  // fitT fitA and fitF
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].windows*(data[i].fitorder+1);
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=sizeof(CMPTYPE);
    size += h_size[i];
  }
    h_size[i] = data[i].windows*(data[i].fitorder+1)*2*sizeof(CMPTYPE);
    size += h_size[i];
  cudaMalloc((void**)&fitT, size);
  cudaMalloc((void**)&fitA, size);
  cudaMemcpy(offsets+PFITS*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);

  CMPTYPE *h_fitT;
  CMPTYPE *h_fitA; 
  unsigned ic, iw;
  for(i=0;i<numIso;i++){
    size = h_size[i]; 
    h_fitT = (CMPTYPE*)malloc(size);
    h_fitA = (CMPTYPE*)malloc(size);
    for(ic=0;ic<=data[i].fitorder;ic++){
      for(iw=0;iw<data[i].windows;iw++){
        h_fitT[ic*data[i].windows+iw] = data[i].fit[findex(iw,ic,FIT_T,data[i].fitorder+1,FIT_F+data[i].fissionable)]; 
     }
    }
    for(ic=0;ic<=data[i].fitorder;ic++){
      for(iw=0;iw<data[i].windows;iw++){
        h_fitA[ic*data[i].windows+iw] = data[i].fit[findex(iw,ic,FIT_A,data[i].fitorder+1,FIT_F+data[i].fissionable)]; 
     }
    }
    gpuErrchk(cudaMemcpy(fitT+h_offset[i],h_fitT,size,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(fitA+h_offset[i],h_fitA,size,cudaMemcpyHostToDevice));
    free(h_fitT);
    free(h_fitA);
  }
//fitF 
  size = 0;
  h_offset[0] = 0;
  for(i=0;i<numIso-1;i++){
    h_size[i] = data[i].fissionable? data[i].windows*(data[i].fitorder+1):0;
    h_offset[i+1]=h_offset[i] + h_size[i];
    h_size[i]*=sizeof(CMPTYPE);
    size += h_size[i];
  }
    h_size[i] = data[i].fissionable? data[i].windows*(data[i].fitorder+1)*sizeof(CMPTYPE):0;
    size += h_size[i];
  cudaMalloc((void**)&fitF, size);
  cudaMemcpy(offsets+PFITF*numIso, h_offset, sizeof(int)*numIso, cudaMemcpyHostToDevice);


  CMPTYPE *h_fitF ;
  for(i=0;i<numIso;i++){
    size = h_size[i]; 
    if(0!=size){
    h_fitF = (CMPTYPE*)malloc(size);
    for(ic=0;ic<=data[i].fitorder;ic++){
      for(iw=0;iw<data[i].windows;iw++){
        h_fitF[ic*data[i].windows+iw] = data[i].fit[findex(iw,ic,FIT_F,data[i].fitorder+1,FIT_F+data[i].fissionable)]; 
     }
    }
    gpuErrchk(cudaMemcpy(fitF+h_offset[i],h_fitF,size,cudaMemcpyHostToDevice));
    free(h_fitF);
    }
  }


  free(h_offset); 
  free(h_size);
 
#if defined(__QUICKWG)
  mtable = wtable;  
#endif
}


multipole::~multipole(){
}

void multipole::release_pointer(){
  gpuErrchk(cudaFree(offsets));
  gpuErrchk(cudaFree(dev_numIso));
  gpuErrchk(cudaFree(dev_integers));
  gpuErrchk(cudaFree(dev_doubles));
  gpuErrchk(cudaFree(mpdata));
  gpuErrchk(cudaFree(l_value));
  gpuErrchk(cudaFree(pseudo_rho));
  gpuErrchk(cudaFree(w_start));
  gpuErrchk(cudaFree(w_end));
  gpuErrchk(cudaFree(fitT));
  gpuErrchk(cudaFree(fitA));
  gpuErrchk(cudaFree(fitF));

#if defined(__QUICKWT)
  unbindwtable();
#endif
}

/*
__device__ void broaden_n_polynomials(double En, double DOPP, double* factors, unsigned n){
//!translated from mit-crpg/WHOPPER
  double sqrtE = sqrt(En);
  double beta  = sqrtE*DOPP;  
  double halfinvDOPP2 = 0.5/(DOPP*DOPP);
  double quarterinvDOPP4 = 0.25/(DOPP*DOPP*DOPP*DOPP);
  double erfBeta, exp_m_beta2;
  if(beta>6.){
    erfBeta = 1.0;
    exp_m_beta2 = 0.0;
  }
  else{
    erfBeta = erf(beta); 
    exp_m_beta2 = exp(-beta*beta);
  }  
  factors[0] = erfBeta/En;
  factors[1] = 1.0/sqrtE;
  factors[2] = erfBeta/En*(halfinvDOPP2+En)+exp_m_beta2/(beta*sqrt(PI));
  if(n>=4){
    factors[3] = 1.0/sqrtE*(En+3.0*halfinvDOPP2);
    for(int i=1;i<n-3;i++){
      factors[i+3] = -factors[i-1]*i*(i+1)*quarterinvDOPP4
                     +factors[i+1]*(En+(3+2*i)*halfinvDOPP2); 
    } 
  }
}
*/
int multipole::findex(int iW, int iC, int type, int orders, int types){
  return iW*orders*types + iC*types + type; 
}

__host__ __device__ int multipole::pindex(int iP, int type){
  return iP*4 + type;
}


