#include "multipole.h"
__host__ __device__  multipole::multipole(struct multipoledata data){
  /*
    TODO:overload = deserves trying; 
    currently, hdf5 is read into multipoledata struct and then
    initialized to multipole object here.
    It is impossible to read hdf5 file directly into this class to be executed on GPU. Maybe make a multipoledata struct in multipole class and overloading = would be faster
  */
  int i;
  fissionable = data.fissionable;
  mode        = data.mode;
  windows     = data.windows;
  sqrtAWR     = data.sqrtAWR;
  startE      = data.startE;
  endE        = data.endE;
  spacing     = data.spacing;
  fitorder    = data.fitorder;
  length      = data.length;
  numL        = data.numL;
  /*  pseudo_rho  = (double*)malloc(numL*sizeof(double));
  for(i=0;i<numL;i++)
    pseudo_rho[i] = data.pseudo_rho[i];
  */
}


__device__  void multipole::xs_eval_fast(double E, double sqrtKT, 
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power, DOPP, DOPP_ECOEF;
  CComplex w_val;

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
  free(twophi);
}

__device__  void multipole::xs_eval_fast(double E,  
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power;
  CComplex PSIIKI, CDUM1, w_val;

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
}


__host__ __device__  int multipole::findex(int type, int iC, int iW){
  return windows*(fitorder+1)*type+windows*iC+iW;
}

__host__ __device__  int multipole::pindex(int type, int iP){
  return length*type + iP;
}

__host__ __device__  void multipole::fill_factors(double sqrtE, double *twophi){
  int iL;
  double arg;
  for(iL = 0; iL<numL; iL++){
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
