#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double
#define CMPTYPE2 double2
#endif

#include "CComplex.h"
#include "multipole.h"
#include <optix.h>

#if defined(__MITW)
#include "Faddeeva.h"
#define w_function Faddeeva::w
#endif

#if defined(__QUICKW)
#include "QuickW.h"
rtBuffer<CComplex<double>, 1>   wtable_buffer;
#endif

#if defined(__FOURIERW)
#include "fourierw.h"
#endif


rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<unsigned int, 1>       dev_integers;
rtBuffer<unsigned int, 1>       dev_numIso;
rtBuffer<unsigned int, 1>       offsets;
rtBuffer<unsigned int, 1>       l_value;

rtBuffer<unsigned int, 1>       w_start;
rtBuffer<unsigned int, 1>       w_end;

rtBuffer<double, 1>       dev_doubles;
rtBuffer<double, 1>       pseudo_rho;

rtBuffer<double, 1>       fitT;
rtBuffer<double, 1>       fitA;
rtBuffer<double, 1>       fitF;

rtBuffer<CComplex<double>, 1>   mpdata;

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

__device__ int pindex(int iP, int type){
  return iP*4 + type;
}


__device__ void fill_factors(int prhoOffset, CMPTYPE sqrtE, int numL,  
                                        CComplex<double> *sigT_factor){
//!translated from mit-crpg/WHOPPER
  int iL;
  double arg;
  double twophi; 
  
  for(iL = 0; iL<numL; iL++){
    twophi = pseudo_rho[iL+prhoOffset] * sqrtE; 
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

RT_CALLABLE_PROGRAM void xs_eval(int iM, CMPTYPE E, CMPTYPE sqrtKT, 
			                 CMPTYPE *sigT, CMPTYPE *sigA, CMPTYPE *sigF){
//!translated from mit-crpg/WHOPPER
  // Copy variables to local memory for efficiency 
  int numIso = dev_numIso[0];
  int tempOffset=iM*DEVREALS;
  // Currently neutrons are slown down from 20.0MeV to 1.0E-5 eV, which is wider than
  // [startE, endE], 
  CMPTYPE startE  = dev_doubles[tempOffset+STARTE];
  CMPTYPE endE    = dev_doubles[tempOffset+ENDE];
  E = (E<startE)*startE + (E>endE)*endE + ((E>=startE)&&(E<=endE))*E;
  CMPTYPE sqrtE = sqrt(E);
  CMPTYPE spacing = dev_doubles[tempOffset+SPACING]; 
  CMPTYPE sqrtAWR = dev_doubles[tempOffset+SQRTAWR];  
  CMPTYPE power, DOPP, DOPP_ECOEF;
  tempOffset = iM*DEVINTS;
  unsigned mode        = dev_integers[tempOffset+MODE];
  unsigned fitorder    = dev_integers[tempOffset+FITORDER];
  unsigned numL        = dev_integers[tempOffset+NUML];
  unsigned fissionable = dev_integers[tempOffset+FISSIONABLE];
  unsigned windows     = dev_integers[tempOffset+WINDOWS];      

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

  tempOffset = offsets[PWIND*numIso+iM]+iW;
  startW = w_start[tempOffset];
  endW   = w_end[tempOffset];
  CComplex<double> sigT_factor[4];
  //CComplex sigtfactor;
  tempOffset = offsets[PPRHO*numIso+iM];
  if(startW <= endW)
    fill_factors(tempOffset,sqrtE,numL,sigT_factor);
  *sigT = 0.0;
  *sigA = 0.0;
  *sigF = 0.0;
  //polynomial fitting

  tempOffset = offsets[PFITS*numIso+iM];
  mode = offsets[PFITF*numIso+iM]; //to save number of registers, use mode as FITF offset since mode has finished life before
  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    *sigT += fitT[tempOffset+iC*windows+iW]*power;
    *sigA += fitA[tempOffset+iC*windows+iW]*power;
    if(MP_FISS == fissionable)
      *sigF += fitF[mode+iC*windows+iW]*power;
 }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

  tempOffset = offsets[PMPDATA*numIso+iM];
  mode       = offsets[PLVAL*numIso+iM];
  for(iP=startW;iP<=endW;iP++){
    //w_val = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP*DOPP_ECOEF;

#if defined(__QUICKWG) 
    w_val =  w_function((sqrtE - mpdata[tempOffset+pindex(iP-1,MP_EA)])*DOPP,mtable)*DOPP_ECOEF;
#else
    w_val =  w_function((sqrtE - mpdata[tempOffset+pindex(iP-1,MP_EA)])*DOPP       )*DOPP_ECOEF;
#endif //end W method


    *sigT += real(mpdata[tempOffset+pindex(iP-1,MP_RT)]*sigT_factor[l_value[mode+iP-1]-1]*w_val);//sigtfactor);	    
    *sigA += real(mpdata[tempOffset+pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      *sigF += real(mpdata[tempOffset+pindex(iP-1,MP_RF)]*w_val);
  }
}


