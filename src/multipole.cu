#include "multipole.h"

#if defined(__QUICKW)
#if defined(__QUICKWG)
// __QUICKWG assigns global memory wtable to multipole member
multipole::multipole(struct multipoledata data, CComplex<CMPTYPE>* wtable){
#else 
// __QUICKWT has bound global memory wtable to texture 
// __QUICKWC declares constant memory as extern in QuickW.cu
multipole::multipole(struct multipoledata data){
#endif //endif __QUICKWG
#else 
// __MITW    uses no wtable
multipole::multipole(struct multipoledata data){
#endif
  size_t size;
  /*
    allocate and assign integers
  */
  size = sizeof(unsigned);
  cudaMalloc((void**)&dev_integers, 5*size);
  cudaMemcpy(dev_integers+MODE,    &(data.mode), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FITORDER, &(data.fitorder), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+NUML, &(data.numL), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FISSIONABLE, &(data.fissionable), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+WINDOWS, &(data.windows),size, cudaMemcpyHostToDevice);

  /*
    allocate and assign doubles
  */
  size = sizeof(CMPTYPE);
  cudaMalloc((void**)&dev_doubles,  3*size);
  cudaMemcpy(dev_doubles+STARTE, &(data.startE), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SPACING,&(data.spacing), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SQRTAWR, &(data.sqrtAWR), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign arrays
  */
  size = data.length*2*sizeof(CMPTYPE);
  CPUComplex<CMPTYPE> *h_ea = (CPUComplex<CMPTYPE>*)malloc(size);
  CPUComplex<CMPTYPE> *h_rt = (CPUComplex<CMPTYPE>*)malloc(size);
  CPUComplex<CMPTYPE> *h_ra = (CPUComplex<CMPTYPE>*)malloc(size);
  CPUComplex<CMPTYPE> *h_rf;
  if(data.fissionable){
  h_rf = (CPUComplex<CMPTYPE>*)malloc(size);
  }

  unsigned ip;
  for(ip=0;ip<data.length;ip++){
    h_ea[ip] = data.mpdata[pindex(ip,MP_EA,MP_RF+data.fissionable)];
    h_rt[ip] = data.mpdata[pindex(ip,MP_RT,MP_RF+data.fissionable)];
    h_ra[ip] = data.mpdata[pindex(ip,MP_RA,MP_RF+data.fissionable)];
    if(data.fissionable){
    h_rf[ip] = data.mpdata[pindex(ip,MP_RF,MP_RF+data.fissionable)];
    }
  }

  cudaMalloc((void**)&mpdata_ea, size);
  cudaMemcpy(mpdata_ea, h_ea, size, cudaMemcpyHostToDevice);
  free(h_ea);

  cudaMalloc((void**)&mpdata_rt, size);
  cudaMemcpy(mpdata_rt, h_rt, size, cudaMemcpyHostToDevice);
  free(h_rt);

  cudaMalloc((void**)&mpdata_ra, size);
  cudaMemcpy(mpdata_ra, h_ra, size, cudaMemcpyHostToDevice);
  free(h_ra);

  if(data.fissionable){
  cudaMalloc((void**)&mpdata_rf, size);
  cudaMemcpy(mpdata_rf, h_rf, size, cudaMemcpyHostToDevice);
  free(h_rf);
  }

  size = data.length*sizeof(unsigned);
  cudaMalloc((void**)&l_value, size);
  cudaMemcpy(l_value, data.l_value, size, cudaMemcpyHostToDevice);

  size = data.numL*sizeof(CMPTYPE);
  cudaMalloc((void**)&pseudo_rho, size);
  cudaMemcpy(pseudo_rho, data.pseudo_rho, size, cudaMemcpyHostToDevice);


  size = data.windows*sizeof(int);
  cudaMalloc((void**)&w_start, size);
  cudaMemcpy(w_start, data.w_start, size, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&w_end, size);
  cudaMemcpy(w_end, data.w_end, size, cudaMemcpyHostToDevice);

  size = (data.fitorder+1)*data.windows*sizeof(CMPTYPE);
  //cudaMalloc((void**)&fit, size);
  //cudaMemcpy(fit, data.fit, size, cudaMemcpyHostToDevice);
  unsigned ic, iw;
  CMPTYPE *h_fitT = (CMPTYPE*)malloc(size);
  CMPTYPE *h_fitA = (CMPTYPE*)malloc(size);
  CMPTYPE *h_fitF; 
  if(data.fissionable){
  h_fitF = (CMPTYPE*)malloc(size);
  }
  for(ic=0;ic<=data.fitorder;ic++){
    for(iw=0;iw<data.windows;iw++){
      h_fitT[ic*data.windows+iw] = data.fit[findex(iw,ic,FIT_T,data.fitorder+1,2+data.fissionable)]; 
   }
  }
  for(ic=0;ic<=data.fitorder;ic++){
    for(iw=0;iw<data.windows;iw++){
      h_fitA[ic*data.windows+iw] = data.fit[findex(iw,ic,FIT_A,data.fitorder+1,2+data.fissionable)]; 
   }
  }
  if(data.fissionable){
  for(ic=0;ic<=data.fitorder;ic++){
    for(iw=0;iw<data.windows;iw++){
      h_fitF[ic*data.windows+iw] = data.fit[findex(iw,ic,FIT_F,data.fitorder+1,2+data.fissionable)]; 
   }
  }
  }

  gpuErrchk(cudaMalloc((void**)&fitT, size));
  gpuErrchk(cudaMemcpy(fitT,h_fitT,size,cudaMemcpyHostToDevice));
  free(h_fitT);
  gpuErrchk(cudaMalloc((void**)&fitA, size));
  gpuErrchk(cudaMemcpy(fitA,h_fitA,size,cudaMemcpyHostToDevice));
  free(h_fitA);
  if(data.fissionable) {
  gpuErrchk(cudaMalloc((void**)&fitF, size));
  gpuErrchk(cudaMemcpy(fitF,h_fitF,size,cudaMemcpyHostToDevice));
  free(h_fitF);
  }
 
 
#if defined(__QUICKWG)
  mtable = wtable;  
#endif
}


multipole::~multipole(){
}

void multipole::release_pointer(unsigned fissionable){
  gpuErrchk(cudaFree(dev_integers));
  gpuErrchk(cudaFree(dev_doubles));
  gpuErrchk(cudaFree(mpdata_ea));
  gpuErrchk(cudaFree(mpdata_rt));
  gpuErrchk(cudaFree(mpdata_ra));
  gpuErrchk(cudaFree(l_value));
  gpuErrchk(cudaFree(pseudo_rho));
  gpuErrchk(cudaFree(w_start));
  gpuErrchk(cudaFree(w_end));
  gpuErrchk(cudaFree(fitT));
  gpuErrchk(cudaFree(fitA));
if(fissionable){
  gpuErrchk(cudaFree(fitF));
  gpuErrchk(cudaFree(mpdata_rf));
}
#if defined(__QUICKWT)
  unbindwtable();
#endif
}

// xs eval with MIT Faddeeva()
#if defined(__MITW) || defined(__QUICKW) || defined(__FOURIERW)
__device__  void multipole::xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, 
			                 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){

  // Copy variables to local memory for efficiency 
  CMPTYPE sqrtE = sqrt(E);
  CMPTYPE spacing = dev_doubles[SPACING];
  CMPTYPE startE  = dev_doubles[STARTE];
  CMPTYPE sqrtAWR = dev_doubles[SQRTAWR];
  CMPTYPE power, DOPP, DOPP_ECOEF;
  unsigned mode        = dev_integers[MODE];
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned numL        = dev_integers[NUML];
  unsigned fissionable = dev_integers[FISSIONABLE];
  unsigned windows     = dev_integers[WINDOWS];

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

  startW = w_start[iW];
  endW   = w_end[iW];
  CComplex<double> sigT_factor[4];
  //CComplex sigtfactor;
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += fitT[iC*windows+iW]*power;
    sigA += fitA[iC*windows+iW]*power;
    if(MP_FISS == fissionable)
      sigF += fitF[iC*windows+iW]*power;
 }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

#if defined(__TRACK)
  numL = 0;
#endif
  for(iP=startW;iP<=endW;iP++){
    //w_val = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP*DOPP_ECOEF;

#if defined(__CFLOAT)
    CComplex<float>  zfloat  = mpdata_ea[iP-1];//mpdata[pindex(iP-1,MP_EA)];
    CComplex<double> zdouble = CComplex<double>((double)real(zfloat),(double)imag(zfloat));

#if defined(__QUICKWG) 
    w_val =  w_function(((double)sqrtE - zdouble)*(double)DOPP,mtable)*(double)DOPP_ECOEF;
#else
    w_val =  w_function(((double)sqrtE - zdouble)*(double)DOPP       )*(double)DOPP_ECOEF;
#endif //end W method

#else //not defined __CFLOAT

#if defined(__QUICKWG) 
    w_val =  w_function((sqrtE - mpdata_ea[iP-1])*DOPP,mtable)*DOPP_ECOEF;
#else
    w_val =  w_function((sqrtE - mpdata_ea[iP-1])*DOPP       )*DOPP_ECOEF;
#endif //end W method

#endif //end if __CFLOAT

#if defined(__TRACK)
    numL++;
#endif 
#if defined(__PLOT)
  if(threadIdx.x<=50){
    CComplex<CMPTYPE> zout = (sqrtE - mpdata_ea[iP-1])*DOPP;
    printf("%+20.16e %+20.16e\n", real(zout),imag(zout));
}
#endif
#if defined(__CFLOAT)
    zfloat = mpdata_rt[iP-1]; 
    zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat))*sigT_factor[l_value[iP-1]-1];
    sigT += (CMPTYPE)real(zdouble*w_val);

    zfloat = mpdata_ra[iP-1]; 
    zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat));
    sigA += (CMPTYPE)real(zdouble*w_val);
    if(MP_FISS == fissionable){
      zfloat = mpdata_rf[iP-1]; 
      zdouble= CComplex<double>((double)real(zfloat),(double)imag(zfloat));
      sigF += (CMPTYPE)real(zdouble*w_val);
    }

#else
    sigT += real(mpdata_rt[iP-1]*sigT_factor[l_value[iP-1]-1]*w_val);//sigtfactor);	    
    sigA += real(mpdata_ra[iP-1]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata_rf[iP-1]*w_val);
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

int multipole::pindex(int iP, int type, int types){
  return iP*types + type;
}

__device__ void multipole::fill_factors(CMPTYPE sqrtE, int numL, 
                                        CComplex<double> *sigT_factor){
  int iL;
  double arg;
  double twophi; 
  
  for(iL = 0; iL<numL; iL++){
    twophi = pseudo_rho[iL] * sqrtE; 
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
