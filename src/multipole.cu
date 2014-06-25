#include "multipole.h"
#if defined(__QUICKWT)
#if defined(__CFLOAT)
texture<float2, 2> tex_wtable;
#else
texture<double2, 2> tex_wtable;
#endif
#endif

#if defined(__QUICKWG) || defined(__QUICKWT)
multipole::multipole(struct multipoledata data, CComplex<CMPTYPE>* wtable){
#else
multipole::multipole(struct multipoledata data){
#endif
  size_t size;
  /*
    allocate and assign integers
  */
  size = sizeof(unsigned);
  cudaMalloc((void**)&dev_integers, 4*size);
  cudaMemcpy(dev_integers+MODE,    &(data.mode), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FITORDER, &(data.fitorder), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+NUML, &(data.numL), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FISSIONABLE, &(data.fissionable), size, cudaMemcpyHostToDevice);

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
  size = data.length*(MP_RF+data.fissionable)*2*sizeof(CMPTYPE);
  cudaMalloc((void**)&mpdata, size);
  cudaMemcpy(mpdata, data.mpdata, size, cudaMemcpyHostToDevice);

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

  size = (FIT_F+data.fissionable)*(data.fitorder+1)*data.windows*sizeof(CMPTYPE);
  cudaMalloc((void**)&fit, size);
  cudaMemcpy(fit, data.fit, size, cudaMemcpyHostToDevice);
  
#if defined(__QUICKWT)
  //cudaBindTexture(NULL, tex_wtable, wtable, LENGTH*LENGTH*sizeof(CMPTYPE)*2);
#if defined(__CFLOAT)
  cudaChannelFormatDesc desc = cudaCreateChannelDesc<float2>();
#else
  cudaChannelFormatDesc desc = cudaCreateChannelDesc<double2>();
#endif
  cudaBindTexture2D(NULL, tex_wtable, wtable, desc, LENGTH, LENGTH, sizeof(CMPTYPE)*2*LENGTH);
#endif

#if defined(__QUICKWG)
  table = wtable;  
#endif

}


multipole::~multipole(){
}

void multipole::release_pointer(){
  gpuErrchk(cudaFree(dev_integers));
  gpuErrchk(cudaFree(dev_doubles));
  gpuErrchk(cudaFree(mpdata));
  gpuErrchk(cudaFree(l_value));
  gpuErrchk(cudaFree(pseudo_rho));
  gpuErrchk(cudaFree(w_start));
  gpuErrchk(cudaFree(w_end));
  gpuErrchk(cudaFree(fit));
#if defined(__QUICKWT)
  cudaUnbindTexture(tex_wtable);
#endif
}

// xs eval with MIT Faddeeva()
#if defined(__MITW) || defined(__QUICKW)
__device__  void multipole::xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, 
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
		       
#if defined(__QUICKWG) 
    w_val = w_function((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP,table)*DOPP_ECOEF;
#endif

#if defined(__QUICKWC) || defined(__QUICKWT)
    w_val = w_function((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP      )*DOPP_ECOEF;
    // __QUICKWT extern texture in QuickW.cu
    // __QUICKWC extern array   in QuickW.cu
#endif


#if defined(__MITW)
#if defined(__CFLOAT)
    CComplex<float>  zfloat  = mpdata[pindex(iP-1,MP_EA)];
    CComplex<double> zdouble = CComplex<double>((double)real(zfloat),(double)imag(zfloat));
    /*CComplex<double> zdouble = CComplex<double>((double)real(mpdata[pindex(iP-1,MP_EA)]),
						(double)imag(mpdata[pindex(iP-1,MP_EA)]));*/
    zdouble = Faddeeva::w(((double)sqrtE - zdouble)*(double)DOPP);
    w_val = CComplex<float>((float)real(zdouble), (float)imag(zdouble))*DOPP_ECOEF;
#else
    w_val = Faddeeva::w((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP,0.0)*DOPP_ECOEF;
#endif
#endif
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*sigT_factor[l_value[iP-1]-1]*w_val);//sigtfactor);	    
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*w_val);
  }

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

__host__ __device__ int multipole::findex(int iW, int iC, int type, int orders, int types){
  return iW*orders*types + iC*types + type; 
}

__host__ __device__ int multipole::pindex(int iP, int type){
  return iP*4 + type;
}

__device__ void multipole::fill_factors(CMPTYPE sqrtE, int numL, 
                                        CComplex<CMPTYPE> *sigT_factor){
  int iL;
  CMPTYPE arg;
  CMPTYPE twophi; 
  
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
    sigT_factor[iL] = CComplex<CMPTYPE>(cos(twophi), -sin(twophi));
  }

}
