#include "multipole.h"
multipole::multipole(struct multipoledata data){
  size_t size;
  /*
    allocate and assign integers
  */
  size = sizeof(int);
  cudaMalloc((void**)&dev_integers, 6*size);
  cudaMemcpy(dev_integers+MODE,    &(data.mode), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+WINDOWS, &(data.windows), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FITORDER, &(data.fitorder), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+NUML, &(data.numL), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+FISSION, &(data.fissionable), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_integers+LENGTH, &(data.length), size, cudaMemcpyHostToDevice);

  /*
    allocate and assign doubles
  */
  size = sizeof(double);
  cudaMalloc((void**)&dev_doubles,  4*size);
  cudaMemcpy(dev_doubles+STARTE, &(data.startE), size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_doubles+ENDE,   &(data.endE), size, cudaMemcpyHostToDevice);
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
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power, DOPP, DOPP_ECOEF;
  CComplex w_val;
  /*
  if(1==mode)
    iW = (int)(sqrtE - sqrt(startE))/spacing;
  else if(2==mode)
    iW = (int)(log(E) - log(startE))/spacing;
  else
    iW = (int)( E - startE )/spacing;
  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE,twophi);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC);
    sigT += fit[findex(FIT_T, iC, iW)]*power;
    sigA += fit[findex(FIT_A, iC, iW)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(FIT_F, iC, iW)]*power;
  }
  //Faddeeva evaluation in advance
  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/sqrt(PI);
  for(iP=startW;iP<=endW;iP++){
    Z_array[iP-startW] = (sqrtE - mpdata[pindex(MP_EA,iP)])*DOPP;
    W_array[iP-startW] = Faddeeva::w(Z_array[iP-startW])*DOPP_ECOEF;
  }

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    sigT += real(mpdata[pindex(MP_RT,iP)]*sigT_factor[l_value[iP]-1]*W_array[iP-startW]);
    sigA += real(mpdata[pindex(MP_RA,iP)]*W_array[iP-startW]);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(MP_RF,iP)]*W_array[iP-startW]);
  }
  */

}

__device__  void multipole::xs_eval_fast(double E,  
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power;
  CComplex PSIIKI, CDUM1, w_val;
  /*
  twophi = (double*)malloc(sizeof(double)*numL);
 
  if(1==mode)
    iW = (int)(sqrtE - sqrt(startE))/spacing;
  else if(2==mode)
    iW = (int)(log(E) - log(startE))/spacing;
  else
    iW = (int)( E - startE )/spacing;
  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE,twophi);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC);
    sigT += fit[findex(FIT_T, iC, iW)]*power;
    sigA += fit[findex(FIT_A, iC, iW)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(FIT_F, iC, iW)]*power;
  }
  //Faddeeva evaluation in advance

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    PSIIKI = -ONEI/(mpdata[pindex(MP_EA,iP)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(MP_RT,iP)]*CDUM1*sigT_factor[l_value[iP]-1]);
    sigA += real(mpdata[pindex(MP_RA,iP)]*CDUM1);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(MP_RF,iP)]*CDUM1);
  }
  free(twophi);
  */
}


__host__ __device__  int multipole::findex(int type, int iC, int iW){
  //  return windows*(fitorder+1)*type+windows*iC+iW;
  return 0;
}

__host__ __device__  int multipole::pindex(int type, int iP){
  //  return length*type + iP;
  return 0;
}

//TODO: here just continue the initilization scheme, it deserves trying make some values shared
__device__ void multipole::fill_factors(double sqrtE, double *twophi, CComplex *sigT_factor){
  int iL;
  double arg;

  for(iL = 0; iL<dev_integers[NUML]; iL++){
    twophi[iL] = pseudo_rho[iL] * sqrtE; 
    if(2==iL)
      twophi[iL] -= atan(twophi[iL]);
    else if(3==iL){
      arg = 3.0*twophi[iL] / (3.0 - twophi[iL]*twophi[iL]);
      twophi[iL] -= atan(arg);
    }
    else if(4==iL){
      arg = twophi[iL]*(15.0 - twophi[iL]*twophi[iL])/(15.0 - 6.0*twophi[iL]*twophi[iL]);
      twophi[iL] -= atan(arg);
    }
    twophi[iL] *= 2.0;
    sigT_factor[iL] = CComplex(cos(twophi[iL]), -sin(twophi[iL]));
  }

}
