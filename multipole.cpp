#include "multipole.h"
multipole::multipole(){

}

multipole::multipole(char filename[]){
  

}

void multipole::xs_eval_fast(double E, double sqrtKT, 
			double *sigT, double *sigA, double *sigF){
  int    iP, iC, iW, startW, endW;
  double *twophi;
  double sqrtE = sqrt(E);
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
