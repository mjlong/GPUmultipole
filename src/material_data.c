#include "material_data.h"

void matread(struct matdata *pmat, char* filename){
  int i,j;
  int m=0;
  unsigned iS;
  float rho;
  FILE *file;
  char line[MAXLEN];
  clearline(line);
  file=fopen(filename,"r");
  //
  //determine size and number of materials
  //
  i=0;
  while(fgets(line,MAXLEN,file)){
    i++;
    for(j=0;j<MAXLEN-1;j++){
      if((' '==line[j]) && (' '!=line[j+1]))
        break;
    }
    if(MAXLEN-1==j)
      m++;
    clearline(line);
  }
  //mat wise data
  pmat->numMat=m;
  pmat->offsets = (unsigned*)malloc(sizeof(unsigned)*m);
  pmat->N_tot = (float*)malloc(sizeof(float)*m);
  printf("mat=%d, niso=%d\n",m,i-m);
  //isotope wise data
  pmat->isotopes = (unsigned*)malloc(sizeof(unsigned)*(i-m));
  pmat->densities= (float*)malloc(sizeof(float)*(i-m));
  rewind(file);
  clearline(line);
  i=0;
  m=0;
  while(fgets(line,MAXLEN,file)){
    for(j=0;j<MAXLEN-1;j++){
      if((' '==line[j]) && (' '!=line[j+1]))
        break;
    }
    if(j<MAXLEN-1){
      //line of iS and rho
      pmat->isotopes[i]  = atoi(line);
      pmat->densities[i] = atof(line+j);
      //printf("i=%d: iS=%d, rho=%f\n",i,mat.isotopes[i],mat.densities[i]);
      i++;
    }
    else{
      //line of iM
      pmat->N_tot[m]  = atof(line);
      pmat->offsets[m++]=i;
      //printf("m=%d,offset=%d\n",m-1,i);
    }
    clearline(line);
  }

  fclose(file);

  return;
}


void clearline(char *line){
  int j;
  for(j=0;j<MAXLEN;j++)
    line[j]=' ';
}
