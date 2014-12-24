#include "multipole_data.h"

extern void h5read(struct multipoledata & pole, char filename[]);
unsigned count_isotopes(char* inputname){
  unsigned i=0;
  char line[FILENAMELEN];
  FILE* fp=NULL;
  fp = fopen(inputname,"r");
  if (fp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n",inputname);
    exit(1);
  }
  while(fgets(line, FILENAMELEN, fp)) 
    i++;
  fclose(fp);
  return i;
}

void isotope_read(char* input, struct multipoledata* isotopes ){
  int numIso=-1;
  char line[MAXISOTOPES][FILENAMELEN];
  FILE *fp = NULL;
  fp = fopen(input,"r");

  if (fp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n",input);
    exit(1);
  }
  int i;
  while(fgets(line[++numIso], FILENAMELEN, fp) != NULL) {
    i=-1;
    printf("%s", line[numIso]);
    while('\n'!=line[numIso][++i]){}
    line[numIso][i]='\0';
    printf("\n"); 
  }
  fclose(fp);
  for(int i=0;i<numIso;i++)  
    h5read(isotopes[i],line[i]);
}


void freeMultipoleData(int numIsos, struct multipoledata* data){
  for(int i=0;i<numIsos;i++){
    free(data[i].fit);
    free(data[i].mpdata);
    free(data[i].l_value);
    free(data[i].pseudo_rho);
    free(data[i].w_start);
    free(data[i].w_end);
  }
  free(data);
}


#if defined(__ALLCPU)
void host_xs_eval_fast(struct multipoledata iso, CMPTYPE E, CMPTYPE sqrtKT, 
			                 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF){
//!translated from mit-crpg/WHOPPER
  // Currently neutrons are slown down from 20.0MeV to 1.0E-5 eV, which is wider than
  // [startE, endE], 
  CMPTYPE startE  = iso.startE;
  CMPTYPE endE    = iso.endE;
  E = (E<startE)*startE + (E>endE)*endE + ((E>=startE)&&(E<=endE))*E;
  CMPTYPE sqrtE = sqrt(E);
  CMPTYPE spacing = iso.spacing; 
  CMPTYPE sqrtAWR = iso.sqrtAWR;  
  CMPTYPE power, DOPP, DOPP_ECOEF;
  unsigned mode        = iso.mode;
  unsigned fitorder    = iso.fitorder;
  unsigned numL        = iso.numL;
  unsigned fissionable = fissionable;
  unsigned windows     = windows;      

  int    iP, iC, iW, startW, endW;
  if(1==mode)
    iW = (int)((sqrtE - sqrt(startE))/spacing);
  else if(2==mode)
    iW = (int)((log(E) - log(startE))/spacing);
  else
    iW = (int)(( E - startE )/spacing);
  //iW =  (int)(((sqrtE - sqrt(startE))*(1==mode) + (log(E) - log(startE))*(2==mode) +  ( E - startE )*(3==mode))/spacing);
  //CPUComplex<CMPTYPE> w_val;
  CPUComplex<double> w_val;

  startW = iso.w_start[iW];
  endW   = iso.w_end[iW];
  CPUComplex<double> sigT_factor[4];
  //CPUComplex sigtfactor;
  if(startW <= endW)
    fill_factors(sqrtE,numL,iso.pseudo_rho,sigT_factor);
  sigT = 0.0;
  sigA = 0.0;
  sigF = 0.0;
  //polynomial fitting

  for (iC=0;iC<=fitorder;iC++){
    power = (CMPTYPE)pow((double)E,(double)iC*0.5-1.0);
    sigT += iso.fit[findex(iW,iC,FIT_T,fitorder+1,FIT_F+fissionable)]*power;
    sigA += iso.fit[findex(iW,iC,FIT_A,fitorder+1,FIT_F+fissionable)]*power;
    if(MP_FISS == fissionable)
      sigF += iso.fit[findex(iW,iC,FIT_F,fitorder+1,FIT_F+fissionable)]*power;
 }

  DOPP = sqrtAWR/sqrtKT;
  DOPP_ECOEF = DOPP/E*sqrt(PI);

  for(iP=startW;iP<=endW;iP++){
    //w_val = (sqrtE - mpdata[pindex(iP-1,MP_EA)])*DOPP*DOPP_ECOEF;

#if defined(__QUICKWG) 
    w_val =  w_function((sqrtE - iso.mpdata[pindex(iP-1,MP_EA)])*DOPP,mtable)*DOPP_ECOEF;
#else
    w_val =  w_function((sqrtE - iso.mpdata[pindex(iP-1,MP_EA)])*DOPP       )*DOPP_ECOEF;
#endif //end W method

    sigT += real(iso.mpdata[pindex(iP-1,MP_RT)]*sigT_factor[iso.l_value[mode+iP-1]-1]*w_val);//sigtfactor);	    
    sigA += real(iso.mpdata[pindex(iP-1,MP_RA)]*w_val);                              
    if(MP_FISS == fissionable)
      sigF += real(iso.mpdata[pindex(iP-1,MP_RF)]*w_val);
  }
}


void fill_factors(CMPTYPE sqrtE, int numL, CMPTYPE* pseudo_rho,   
                                        CPUComplex<double> *sigT_factor){
//!translated from mit-crpg/WHOPPER
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
    sigT_factor[iL] = CPUComplex<double>(cos(twophi), -sin(twophi));
  }
}

int pindex(int iP, int type){
  return iP*4 + type;
}
int findex(int iW, int iC, int type, int orders, int types){
  return iW*orders*types + iC*types + type; 
}
#endif
