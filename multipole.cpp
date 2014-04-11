#include "multipole.h"
multipole::multipole(){

}

//TODO: this doesn't work
/*
multipole::multipole(char filename[]){
  isotope iso(filename);
  isotopeinfo(iso);

}

void multipole::isotopeinfo(isotope iso){
  atomic_weight_ratio = iso.awr;
  pseudo_rho = (double*)malloc(sizeof(double)*iso.number_l);
  for(int i=0;i<iso.number_l;i++)
    pseudo_rho[i]=iso.pseudo_rho[i];
}
*/

void multipole::xs_eval_fast(double E, double sqrtKT, 
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW, cnt,maxwindow=0;
  double *twophi;
  double sqrtE = sqrt(E);
  double power, DOPP, DOPP_ECOEF;
  complex<double> PSIIKI, CDUM1, w_val;
  complex<double> *Z_array, *W_array;

  twophi = (double*)malloc(sizeof(double)*numL);
  for(iW=0;iW<windows;iW++){
    cnt = w_end[iW]-w_start[iW] + 1;
    if(cnt > maxwindow)
      maxwindow = cnt;
  }
  Z_array = (complex<double>*)malloc(sizeof(complex<double>)*maxwindow);
  W_array = (complex<double>*)malloc(sizeof(complex<double>)*maxwindow);
  
  if(1==mode)
    iW = (int)(sqrtE - sqrt(startE))/spacing;
  else if(2==mode)
    iW = (int)(log(E) - log(startE))/spacing;
  else
    iW = (int)( E - startE )/spacing;
  startW = w_start[iW];
  endW   = w_end[iW];
  if(startW <= endW)
    fill_factors(sqrtE);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<fitorder;iC++){
    power = pow(E,iC);
    sigT += fit[findex(FIT_T, iC, iW)]*power;
    sigA += fit[findex(FIT_A, iC, iW)]*power;
    if(MP_FISS == fissionable)
      sigF += fit[findex(FIT_F, iC, iW)]*power;
  }
  //Faddeeva evaluation in advance
  DOPP = sqrt(atomic_weight_ratio)/sqrtKT;
  DOPP_ECOEF = DOPP/sqrt(PI);
  for(iP=startW;iP<=endW;iP++){
    Z_array[iP-startW] = (sqrtE - mpdata[pindex(MP_EA,iP)])*DOPP;
    W_array[iP-startW] = w(Z_array[iP-startW])*DOPP_ECOEF;
  }

  //evaluating
  for(iP=startW;iP<=endW;iP++){
    sigT += real(mpdata[pindex(MP_RT,iP)]*sigT_factor[l_value[iP]-1]*W_array[iP-startW]);
    sigA += real(mpdata[pindex(MP_RA,iP)]*W_array[iP-startW]);
    if(MP_FISS == fissionable)
      sigF += real(mpdata[pindex(MP_RF,iP)]*W_array[iP-startW]);
  }
}

int multipole::findex(int type, int iC, int iW){
  return windows*(fitorder+1)*type+windows*iC+iW;
}

int multipole::pindex(int type, int iP){
  return length*type + iP;
}

void multipole::fill_factors(double sqrtE){
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
    sigT_factor[iL] = complex<double>(cos(twophi[iL]),-sin(twophi[iL]));
  }

}
