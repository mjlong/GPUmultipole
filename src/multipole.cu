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
    cudaMemcpy(dev_numIso, &numIso, sizeof(int), cudaMemcpyHostToDevice);
  /*
    allocate and assign doubles
  */
  size = sizeof(CMPTYPE);
  cudaMalloc((void**)&dev_doubles,  DEVREALS*size*numIso);
  for(i=0;i<numIso;i++){
    cudaMemcpy(dev_doubles+i*DEVINTS+STARTE,  &(data[i].startE),  size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_doubles+i*DEVINTS+SPACING ,&(data[i].spacing), size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_doubles+i*DEVINTS+SQRTAWR, &(data[i].sqrtAWR), size, cudaMemcpyHostToDevice);
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
    gpuErrchk(cudaMemcpy(fitA+h_offset[i],h_fitT,size,cudaMemcpyHostToDevice));
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

// xs eval with MIT Faddeeva()
#if defined(__MITW) || defined(__QUICKW) || defined(__FOURIERW)
__device__  void multipole::xs_eval_fast(int iM, CMPTYPE E, CMPTYPE sqrtKT, 
			                 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){

  // Copy variables to local memory for efficiency 
  int numIso = dev_numIso[0];
  int tempOffset=iM*numIso;
  CMPTYPE sqrtE = sqrt(E);
  CMPTYPE spacing = dev_doubles[tempOffset+SPACING]; 
  CMPTYPE startE  = dev_doubles[tempOffset+STARTE];
  CMPTYPE sqrtAWR = dev_doubles[tempOffset+SQRTAWR];  
  CMPTYPE power, DOPP, DOPP_ECOEF;
  unsigned mode        = dev_integers[iM*numIso+MODE];
  unsigned fitorder    = dev_integers[iM*numIso+FITORDER];
  unsigned numL        = dev_integers[iM*numIso+NUML];
  unsigned fissionable = dev_integers[iM*numIso+FISSIONABLE];
  unsigned windows     = dev_integers[iM*numIso+WINDOWS];      //so far tempOffset = iM*numIso

  int    iP, iC, iW, startW, endW;
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  //iW =  (int)(((sqrtE - sqrt(startE))*(1==mode) + (log(E) - log(startE))*(2==mode) +  ( E - startE )*(3==mode))/spacing);
  //CComplex<CMPTYPE> w_val;
  CComplex<double> w_val;

  tempOffset = offsets[PWIND*numIso+iM]+iW;
  startW = w_start[tempOffset];
  endW   = w_end[tempOffset];
  CComplex<double> sigT_factor[4];
  //CComplex sigtfactor;
  tempOffset = offsets[PLVAL*numIso+iM];
  if(startW <= endW)
    fill_factors(tempOffset,sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  tempOffset = offsets[PFITS*numIso+iM];
  mode = offsets[PFITF*numIso+iM]; //to save number of registers, use mode as FITF offset since mode has finished life before
  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += fitT[tempOffset+iC*windows+iW]*power;
    sigA += fitA[tempOffset+iC*windows+iW]*power;
    if(MP_FISS == fissionable)
      sigF += fitF[mode+iC*windows+iW]*power;
 }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

#if defined(__TRACK)
  numL = 0;
#endif
  tempOffset = offsets[PMPDATA*numIso+iM];
  mode       = offsets[PLVAL*numIso+iM];
  for(iP=startW;iP<=endW;iP++){
    //w_val = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP*DOPP_ECOEF;

#if defined(__CFLOAT)
    CComplex<float>  zfloat  = mpdata[tempOffset+pindex(iP-1,MP_EA)];
    CComplex<double> zdouble = CComplex<double>((double)real(zfloat),(double)imag(zfloat));

#if defined(__QUICKWG) 
    w_val =  w_function(((double)sqrtE - zdouble)*(double)DOPP,mtable)*(double)DOPP_ECOEF;
#else
    w_val =  w_function(((double)sqrtE - zdouble)*(double)DOPP       )*(double)DOPP_ECOEF;
#endif //end W method

#else //not defined __CFLOAT

#if defined(__QUICKWG) 
    w_val =  w_function((sqrtE - mpdata[tempOffset+pindex(iP-1,MP_EA)])*DOPP,mtable)*DOPP_ECOEF;
#else
    w_val =  w_function((sqrtE - mpdata[tempOffset+pindex(iP-1,MP_EA)])*DOPP       )*DOPP_ECOEF;
#endif //end W method

#endif //end if __CFLOAT

#if defined(__TRACK)
    numL++;
#endif 
#if defined(__PLOT)
  if(threadIdx.x<=50){
    CComplex<CMPTYPE> zout = (sqrtE - mpdata[tempOffset+pindex(iP-1,MP_EA)])*DOPP;
    printf("%+20.16e %+20.16e\n", real(zout),imag(zout));
}
#endif
#if defined(__CFLOAT)
    zfloat = mpdata[pindex(iP-1,MP_RT)]; 
    zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat))*sigT_factor[l_value[mode+iP-1]-1];
    sigT += (CMPTYPE)real(zdouble*w_val);

    zfloat = mpdata[pindex(iP-1,MP_RA)]; 
    zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat));
    sigA += (CMPTYPE)real(zdouble*w_val);
    if(MP_FISS == fissionable){
      zfloat = mpdata[pindex(iP-1,MP_RF)]; 
      zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat));
      sigF += (CMPTYPE)real(zdouble*w_val);
    }

#else
    sigT += real(mpdata[tempOffset+pindex(iP-1,MP_RT)]*sigT_factor[l_value[mode+iP-1]-1]*w_val);//sigtfactor);	    
    sigA += real(mpdata[tempOffset+pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata[tempOffset+pindex(iP-1,MP_RF)]*w_val);
#endif
  }
#if defined(__TRACK)
  sigF = 1.0*numL;
#endif

}
#endif

//xs eval with Quick W()
/*
__device__  void multipole::xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, CComplex<CMPTYPE> *table, 
			                 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){

  // Copy variables to local memory for efficiency 
  unsigned mode        = dev_integers[MODE];
  int    iP, iC, iW, startW, endW;
  CMPTYPE spacing = dev_doubles[SPACING];
  CMPTYPE startE  = dev_doubles[STARTE];
  CMPTYPE sqrtE = sqrt(E);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned numL        = dev_integers[NUML];
  unsigned fissionable = dev_integers[FISSIONABLE];

  CMPTYPE sqrtAWR = dev_doubles[SQRTAWR];
  CMPTYPE power, DOPP, DOPP_ECOEF;
  CComplex<CMPTYPE> w_val;

  startW = w_start[iW];
  endW   = w_end[iW];
  CComplex<CMPTYPE> sigT_factor[4];
  //CComplex sigtfactor;
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,2+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,2+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,2+fissionable)]*power;
  }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

  for(iP=startW;iP<=endW;iP++){
    //sigtfactor = sigT_factor[l_value[iP-1]-1];
    //w_val = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP*DOPP_ECOEF;
    w_val = w_function((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP,table)*DOPP_ECOEF;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*sigT_factor[l_value[iP-1]-1]*w_val);//sigtfactor);	    
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*w_val);
  }

}
*/

//xs eval at 0K
#if defined(__SAMPLE)
__device__  void multipole::xs_eval_fast(CMPTYPE E,  
                        	 	 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){

  // Copy variables to local memory for efficiency 
  unsigned mode        = dev_integers[MODE];
  int    iP, iC, iW, startW, endW;
  CMPTYPE spacing = dev_doubles[SPACING];
  CMPTYPE startE  = dev_doubles[STARTE];
  CMPTYPE sqrtE = sqrt(E);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned fissionable = dev_integers[FISSIONABLE];
  unsigned numL        = dev_integers[NUML];

  CMPTYPE power;
  CComplex<CMPTYPE> PSIIKI, CDUM1, w_val;

 
  startW = w_start[iW];
  endW   = w_end[iW];
  CComplex<CMPTYPE> sigT_factor[4];
  //CComplex sigtfactor;
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,FIT_F+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,FIT_F+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,FIT_F+fissionable)]*power;
  }




  for(iP=startW;iP<=endW;iP++){
    //sigtfactor = sigT_factor[l_value[iP-1]-1];
    PSIIKI = -ONEI/(mpdata[pindex(iP-1,MP_EA)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*CDUM1*sigT_factor[l_value[iP-1]-1]);//sigtfactor);
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*CDUM1);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*CDUM1);
  }
  
}
#endif
//xs eval at 0k but sampled to sqrtKT
/*
__device__  void multipole::xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, CMPTYPE rnd, 
                        	 	 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){

  // Copy variables to local memory for efficiency 
  unsigned mode        = dev_integers[MODE];
  int    iP, iC, iW, startW, endW;
  CMPTYPE spacing = dev_doubles[SPACING];
  CMPTYPE startE  = dev_doubles[STARTE];
  CMPTYPE sqrtAWR = dev_doubles[SQRTAWR];

  E = E + rnd * sqrtKT * sqrt(0.5) / sqrtAWR;
  CMPTYPE sqrtE = sqrt(E);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned fissionable = dev_integers[FISSIONABLE];
  unsigned numL        = dev_integers[NUML];

  CMPTYPE power;
  CComplex<CMPTYPE> PSIIKI, CDUM1, w_val;
 
  startW = w_start[iW];
  endW   = w_end[iW];
  CComplex<CMPTYPE> sigT_factor[4];
  //CComplex sigtfactor;
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,2+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,2+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,2+fissionable)]*power;
  }




  for(iP=startW;iP<=endW;iP++){
    //sigtfactor = sigT_factor[l_value[iP-1]-1];
    PSIIKI = -ONEI/(mpdata[pindex(iP-1,MP_EA)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*CDUM1*sigT_factor[l_value[iP-1]-1]);//sigtfactor);
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*CDUM1);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*CDUM1);
  }
  
}
*/

int multipole::findex(int iW, int iC, int type, int orders, int types){
  return iW*orders*types + iC*types + type; 
}

__host__ __device__ int multipole::pindex(int iP, int type){
  return iP*4 + type;
}

__device__ void multipole::fill_factors(int prhoOffset, CMPTYPE sqrtE, int numL,  
                                        CComplex<double> *sigT_factor){
  int iL;
  double arg;
  double twophi; 
  
  for(iL = 0; iL<numL; iL++){
    twophi = pseudo_rho[iL+prhoOffset] * sqrtE; 
    if(1==iL)
      twophi -= atan(twophi);
    else if(2==iL){
      arg = 3.0*twophi / (3.0 - twophi*twophi);
      twophi -= atan(arg);
    }
    else if(3==iL){
      arg = twophi*(15.0 - twophi*twophi)/(15.0 - 6.0*twophi*twophi);
      twophi -= atan(arg);
    }
    twophi *= 2.0;
    sigT_factor[iL] = CComplex<double>(cos(twophi), -sin(twophi));
  }

}
