#include "multipole_data.h"

extern void h5read(struct multipoledata & pole, char filename[]);

unsigned generateZ(struct multipoledata data, double sqrtKT, double* energy, unsigned num, CPUComplex<CMPTYPE> **pz){
  CMPTYPE E;
  int iP,iC,iW,startW,endW;
  unsigned numz;
  unsigned *sw=(unsigned*)malloc(sizeof(unsigned)*num);
  unsigned *ew=(unsigned*)malloc(sizeof(unsigned)*num);
  
  numz = 0;
  for(int i=0;i<num;i++){
    E = energy[i];
    E = (E<data.startE)*data.startE + (E>data.endE)*data.endE + ((E>=data.startE)&&(E<=data.endE))*E;
 
    if(1==data.mode)
      iW = (int)((sqrt(E) - sqrt(data.startE))/data.spacing);
    else if(2==data.mode)
      iW = (int)((log(E) - log(data.startE))/data.spacing);
    else 
      iW = (int)(( E - data.startE )/data.spacing);
    startW = data.w_start[iW];
    endW = data.w_end[iW];
    numz+=endW-startW+1;
    sw[i] = startW;
    ew[i] = endW;
  }

  *pz = (CPUComplex<CMPTYPE>*)malloc(sizeof(CMPTYPE)*2*numz);
  unsigned iz=0;
  for(int i=0;i<num;i++){
    E = energy[i];
    E = (E<data.startE)*data.startE + (E>data.endE)*data.endE + ((E>=data.startE)&&(E<=data.endE))*E;
    for(iP=sw[i];iP<=ew[i];iP++){
      (*pz)[iz++] = (sqrt(E) - data.mpdata[4*(iP-1)+0])*data.sqrtAWR/sqrtKT;
    }
  }
  free(sw);
  free(ew);
  return numz; 
} 


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





