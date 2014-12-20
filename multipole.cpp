#include "multipole.h"
multipole::multipole(){

}


void multipole::xs_eval_fast(double E, double sqrtKT, 
			     double &sigT, double &sigA, double &sigF, unsigned *counts,unsigned edge){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power, DOPP, DOPP_ECOEF;
  CComplex w_val;

  twophi = (double*)malloc(sizeof(double)*numL);
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);

  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE,twophi);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<=fitorder;iC++){
    power = pow(E,iC*0.5-1.0);
    sigT += fit[findex(iW,iC,FIT_T)]*power;
    sigA += fit[findex(iW,iC,FIT_A)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F)]*power;
  }
  //Faddeeva evaluation in advance
  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);
  int row,col;
  for(iP=startW;iP<=endW;iP++){
    Z_array[iP-startW] = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP;
    CComplex temp = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP;
    col = real(temp)/DELTA+N;
    if(col<0)
      col = 0;
    if(col>=2*N)
      col = 2*N-1;
    row = imag(temp)/DELTA+N;
    if(row<0)
      row = 0;
    if(row>=2*N)
      row = 2*N-1;
    counts[row*edge+col]++;
    W_array[iP-startW] = w(Z_array[iP-startW])*DOPP_ECOEF;
  }

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*sigT_factor[l_value[iP-1]-1]*W_array[iP-startW]);
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*W_array[iP-startW]);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*W_array[iP-startW]);
  }
  free(twophi);
}

void multipole::xs_eval_fast(double E,  
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power;
  CComplex PSIIKI, CDUM1, w_val;

  twophi = (double*)malloc(sizeof(double)*numL);
 
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
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
    sigT += fit[findex(iW,iC,FIT_T)]*power;
    sigA += fit[findex(iW,iC,FIT_A)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(iW,iC,FIT_F)]*power;
  }
  //Faddeeva evaluation in advance

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    PSIIKI = -ONEI/(mpdata[pindex(iP-1,MP_EA)] - sqrtE);
    CDUM1  = PSIIKI / E;
    sigT += real(mpdata[pindex(iP-1,MP_RT)]*CDUM1*sigT_factor[l_value[iP-1]-1]);
    sigA += real(mpdata[pindex(iP-1,MP_RA)]*CDUM1);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(iP-1,MP_RF)]*CDUM1);
  }
  free(twophi);
}


int multipole::findex(int iW, int iC, int type){
  return iW*(fitorder+1)*(2+fissionable) + iC*(2+fissionable) + type;
}

int multipole::pindex(int iP, int type){
  return iP*4 + type;
}

void multipole::fill_factors(double sqrtE, double *twophi){
  int iL;
  double arg;
  for(iL = 0; iL<numL; iL++){
    twophi[iL] = pseudo_rho[iL] * sqrtE; 
    if(1==iL)
      twophi[iL] -= atan(twophi[iL]);
    else if(2==iL){
      arg = 3.0*twophi[iL] / (3.0 - twophi[iL]*twophi[iL]);
      twophi[iL] -= atan(arg);
    }
    else if(3==iL){
      arg = twophi[iL]*(15.0 - twophi[iL]*twophi[iL])/(15.0 - 6.0*twophi[iL]*twophi[iL]);
      twophi[iL] -= atan(arg);
    }
    twophi[iL] *= 2.0;
    sigT_factor[iL] = CComplex(cos(twophi[iL]), -sin(twophi[iL]));
  }

}
