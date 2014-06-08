#include "multipole.h"
texture<unsigned> texl_value;
texture<int2> texfit;
texture<int4> texmpdata;

static __inline__ __device__ double tex1Dfetch_double(texture<int2> t, int i){
  int2 v = tex1Dfetch(t,i);
  return __hiloint2double(v.y, v.x);
}

static __inline__ __device__ CComplex<double> tex1Dfetch_complex(texture<int4> t, int i){
  int4 v = tex1Dfetch(t,i);
  return CComplex<double>(__hiloint2double(v.y, v.x),__hiloint2double(v.w,v.z));
}

multipole::multipole(struct multipoledata data){
  size_t size;
  /*
    allocate and assign integers
  */
  size = sizeof(unsigned);
  gpuErrchk(cudaMalloc((void**)&dev_integers, 4*size));
  cudaMemcpy(dev_integers+MODE,    &(data.mode), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FITORDER, &(data.fitorder), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+NUML, &(data.numL), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FISSIONABLE, &(data.fissionable), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign doubles
  */
  size = sizeof(SETTYPE);
  gpuErrchk(cudaMalloc((void**)&dev_doubles,  3*size));
  cudaMemcpy(dev_doubles+STARTE, &(data.startE), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SPACING,&(data.spacing), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SQRTAWR, &(data.sqrtAWR), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign arrays
  */
  size = data.length*(MP_RF+data.fissionable)*2*sizeof(SETTYPE);
  gpuErrchk(cudaMalloc((void**)&mpdata, size));
  cudaMemcpy(mpdata, data.mpdata, size, cudaMemcpyHostToDevice);
  cudaBindTexture(NULL,texmpdata,mpdata,size);

  size = data.length*sizeof(unsigned);
  gpuErrchk(cudaMalloc((void**)&l_value, size));
  cudaMemcpy(l_value, data.l_value, size, cudaMemcpyHostToDevice);
  cudaBindTexture(NULL,texl_value, l_value, size);

  size = data.numL*sizeof(SETTYPE);
  gpuErrchk(cudaMalloc((void**)&pseudo_rho, size));
  cudaMemcpy(pseudo_rho, data.pseudo_rho, size, cudaMemcpyHostToDevice);


  size = data.windows*sizeof(unsigned);
  gpuErrchk(cudaMalloc((void**)&w_start, size));
  cudaMemcpy(w_start, data.w_start, size, cudaMemcpyHostToDevice);
  gpuErrchk(cudaMalloc((void**)&w_end, size));
  cudaMemcpy(w_end, data.w_end, size, cudaMemcpyHostToDevice);
  //cudaBindTexture(NULL,dtex.W_start, w_start, size);
  //cudaBindTexture(NULL,dtex.W_end,   w_end,   size);

  size = (FIT_F+data.fissionable)*(data.fitorder+1)*data.windows*sizeof(SETTYPE);
  gpuErrchk(cudaMalloc((void**)&fit, size));
  cudaMemcpy(fit, data.fit, size, cudaMemcpyHostToDevice);
  cudaBindTexture(NULL,texfit, fit, size);
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
  //cudaUnbindTexture(dtex.W_start);
  //cudaUnbindTexture(dtex.W_end);
  cudaUnbindTexture(texl_value);
  cudaUnbindTexture(texfit);
  cudaUnbindTexture(texmpdata);
}
__device__  void multipole::xs_eval_fast(SETTYPE E, SETTYPE sqrtKT, 
			                 SETTYPE &sigT, SETTYPE &sigA, SETTYPE &sigF){

  /* Copy variables to local memory for efficiency */ 
  unsigned mode        = dev_integers[MODE];
  int    iP, iC, iW, startW, endW;
  SETTYPE spacing = dev_doubles[SPACING];
  SETTYPE startE  = dev_doubles[STARTE];
  SETTYPE sqrtE = sqrt(E);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned numL        = dev_integers[NUML];
  unsigned fissionable = dev_integers[FISSIONABLE];

  SETTYPE sqrtAWR = dev_doubles[SQRTAWR];
  SETTYPE power, DOPP, DOPP_ECOEF;
  CComplex<SETTYPE> w_val;

  //startW = tex1Dfetch(dtex.W_start,iW);
  //endW   = tex1Dfetch(dtex.W_end,iW);
  startW = w_start[iW];
  endW   = w_end[iW];
  CComplex<SETTYPE> sigT_factor[4];
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC*0.5-1.0);
    //texture
    /*
    sigT += tex1Dfetch_double(texfit,findex(iW,iC,FIT_T,fitorder+1,2+fissionable))*power;
    sigA += tex1Dfetch_double(texfit,findex(iW,iC,FIT_A,fitorder+1,2+fissionable))*power;
    if(MP_FISS == fissionable)
      sigF += tex1Dfetch_double(texfit,findex(iW,iC,FIT_F,fitorder+1,2+fissionable))*power;
    */
    //global
    ///*
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,2+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,2+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,2+fissionable)]*power;
    //*/
  }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

  for(iP=startW;iP<=endW;iP++){
    //texture
    /*
    w_val = Faddeeva::w((sqrtE - tex1Dfetch_complex(texmpdata,pindex(iP-1,MP_EA)))*DOPP)*DOPP_ECOEF;
    sigT += real(tex1Dfetch_complex(texmpdata,pindex(iP-1,MP_RT))*sigT_factor[tex1Dfetch(texl_value,iP-1)-1]*w_val);	    
    sigT += real(tex1Dfetch_complex(texmpdata,pindex(iP-1,MP_RT))*sigT_factor[tex1Dfetch(texl_value,iP-1)-1]*w_val);	    
    sigA += real(tex1Dfetch_complex(texmpdata,pindex(iP-1,MP_RA))*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(tex1Dfetch_complex(texmpdata,pindex(iP-1,MP_RF))*w_val);
    */

    //global
    ///*
    w_val = Faddeeva::w((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP)*DOPP_ECOEF;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*sigT_factor[l_value[iP-1]-1]*w_val);	    
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*w_val);
    //*/
  }

}

__device__  void multipole::xs_eval_fast(SETTYPE E,  
                        	 	 SETTYPE &sigT, SETTYPE &sigA, SETTYPE &sigF){

  /* Copy variables to local memory for efficiency */ 
  unsigned mode        = dev_integers[MODE];
  unsigned fitorder    = dev_integers[FITORDER];
  unsigned fissionable = dev_integers[FISSIONABLE];
  unsigned numL        = dev_integers[NUML];

  SETTYPE spacing = dev_doubles[SPACING];
  SETTYPE startE  = dev_doubles[STARTE];
  
  int    iP, iC, iW, startW, endW;
  CComplex<SETTYPE> sigT_factor[4];
  SETTYPE sqrtE = sqrt(E);
  SETTYPE power;
  CComplex<SETTYPE> PSIIKI, CDUM1, w_val;

  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC);
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,2+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,2+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,2+fissionable)]*power;
  }

  for(iP=startW;iP<=endW;iP++){
    PSIIKI = -ONEI/(mpdata[pindex(iP-1,MP_EA)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*CDUM1*sigT_factor[l_value[iP-1]-1]);
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*CDUM1);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*CDUM1);
  }
  
}



__host__ __device__ int multipole::findex(int iW, int iC, int type, int orders, int types){
  return iW*orders*types + iC*types + type; 
}

__host__ __device__ int multipole::pindex(int iP, int type){
  return iP*4 + type;
}

__device__ void multipole::fill_factors(SETTYPE sqrtE, int numL, 
                                        CComplex<SETTYPE> *sigT_factor){
  int iL;
  SETTYPE arg;
  SETTYPE twophi; 
  
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
    sigT_factor[iL] = CComplex<SETTYPE>(cos(twophi), -sin(twophi));
  }

}
