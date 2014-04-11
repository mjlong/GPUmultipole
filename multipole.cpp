#include "multipole.h"
multipole::multipole(){

}

//TODO: this doesn't work
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

void multipole::xs_eval_fast(double E, double sqrtKT, 
			double &sigT, double &sigA, double &sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
  double power;
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
    fill_factors(sqrtE);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting
  for (iC=0;iC<fitorder;iC++){
    power = pow(E,iC);
    sigT += fit[findex(FIT_T, iC, iW)]*power;
    sigA += fit[findex(FIT_A, iC, iW)]*power;
    if(true == fissionable)
      sigF += fit[findex(FIT_F, iC, iW)]*power;
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
