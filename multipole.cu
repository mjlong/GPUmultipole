#include "multipole.h"
multipole::multipole(struct multipoledata data){
  size_t size;
  /*
    allocate and assign integers
  */
  size = sizeof(int);
  cudaMalloc((void**)&dev_integers, 4*size);
  cudaMemcpy(dev_integers+MODE,    &(data.mode), size, cudaMemcpyHostToDevice);
  //cudaMemcpy(dev_integers+WINDOWS, &(data.windows), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FITORDER, &(data.fitorder), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+NUML, &(data.numL), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FISSIONABLE, &(data.fissionable), size, cudaMemcpyHostToDevice);
  //cudaMemcpy(dev_integers+LENGTH, &(data.length), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign doubles
  */
  size = sizeof(double);
  cudaMalloc((void**)&dev_doubles,  3*size);
  cudaMemcpy(dev_doubles+STARTE, &(data.startE), size, cudaMemcpyHostToDevice);
  //cudaMemcpy(dev_doubles+ENDE,   &(data.endE), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SPACING,&(data.spacing), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+SQRTAWR, &(data.sqrtAWR), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign arrays
  */
  size = data.length*(MP_RF+data.fissionable)*2*sizeof(double);
  cudaMalloc((void**)&mpdata, size);
  cudaMemcpy(mpdata, data.mpdata, size, cudaMemcpyHostToDevice);

  size = data.length*sizeof(unsigned);
  cudaMalloc((void**)&l_value, size);
  cudaMemcpy(l_value, data.l_value, size, cudaMemcpyHostToDevice);

  size = data.numL*sizeof(double);
  cudaMalloc((void**)&pseudo_rho, size);
  cudaMemcpy(pseudo_rho, data.pseudo_rho, size, cudaMemcpyHostToDevice);


  size = data.windows*sizeof(int);
  cudaMalloc((void**)&w_start, size);
  cudaMemcpy(w_start, data.w_start, size, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&w_end, size);
  cudaMemcpy(w_end, data.w_end, size, cudaMemcpyHostToDevice);

  size = (FIT_F+data.fissionable)*(data.fitorder+1)*data.windows*sizeof(double);
  cudaMalloc((void**)&fit, size);
  cudaMemcpy(fit, data.fit, size, cudaMemcpyHostToDevice);
}


multipole::~multipole(){
  cudaFree(dev_integers);
  cudaFree(dev_doubles);
  cudaFree(mpdata);
  cudaFree(l_value);
  cudaFree(pseudo_rho);
  cudaFree(w_start);
  cudaFree(w_end);
  cudaFree(fit);
}
__device__  void multipole::xs_eval_fast(double E, double sqrtKT, 
					 double &sigT, double &sigA, double &sigF, 
					 CComplex* sigT_factor, unsigned blocksize){
  /* Copy variables to local memory for efficiency */ 
  int mode        = dev_integers[MODE];
  int    iP, iC, iW, startW, endW;
  double spacing = dev_doubles[SPACING];
  double startE  = dev_doubles[STARTE];
  double sqrtE = sqrt(E);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  int fitorder    = dev_integers[FITORDER];
  int numL        = dev_integers[NUML];
  int fissionable = dev_integers[FISSIONABLE];
  //int length      = dev_integers[LENGTH];
  //int windows     = dev_integers[WINDOWS];
  //TODO:if length,windows are really not needed, remove them from dev_integers[] array

  //TODO:I've not found wat to allocate for a thread
  // 4 = maximum numL, consistent with max 3==iL in fill_factors()
  //double twophi[4];
  //CComplex sigT_factor[4];

  double sqrtAWR = dev_doubles[SQRTAWR];

  
  double power, DOPP, DOPP_ECOEF;
  CComplex w_val;

  startW = w_start[iW];
  endW   = w_end[iW];

  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor, blocksize);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC*0.5-1.0);
    sigT += fit[findex(iW,iC,FIT_T,fitorder+1,2+fissionable)]*power;
    sigA += fit[findex(iW,iC,FIT_A,fitorder+1,2+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F,fitorder+1,2+fissionable)]*power;
  }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

  for(iP=startW;iP<=endW;iP++){
    w_val = Faddeeva::w((sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP)*DOPP_ECOEF;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*sigT_factor[(l_value[iP-1]-1)*blocksize]*w_val);	    
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*w_val);
  }

}

__device__  void multipole::xs_eval_fast(double E,  
					 double &sigT, double &sigA, double &sigF,
					 CComplex *sigT_factor, unsigned blocksize){
  /* Copy variables to local memory for efficiency */ 
  int mode        = dev_integers[MODE];
  int fitorder    = dev_integers[FITORDER];
  int fissionable = dev_integers[FISSIONABLE];
  //int length      = dev_integers[LENGTH];
  //int windows     = dev_integers[WINDOWS];
  int numL        = dev_integers[NUML];
  //size_t size;
  double spacing = dev_doubles[SPACING];
  double startE  = dev_doubles[STARTE];
  
  int    iP, iC, iW, startW, endW;


  double sqrtE = sqrt(E);
  double power;
  CComplex PSIIKI, CDUM1, w_val;

  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE,numL,sigT_factor, blocksize);
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
  //Faddeeva evaluation in advance

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    PSIIKI = -ONEI/(mpdata[pindex(iP-1,MP_EA)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*CDUM1*sigT_factor[(l_value[iP-1]-1)*blocksize]);
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

//TODO: here just continue the initilization scheme, it deserves trying make some values shared
__device__ void multipole::fill_factors(double sqrtE, int numL, CComplex *sigT_factor, unsigned blocksize){
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
    sigT_factor[iL*blocksize] = CComplex(cos(twophi), -sin(twophi));
  }

}
