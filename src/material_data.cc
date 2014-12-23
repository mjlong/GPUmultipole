#include "material_data.h"

unsigned matread(struct matdata *pmat, char* filename){
  int i,j;
  int m=0;
  unsigned nums;
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
  pmat->offsets = (unsigned*)malloc(sizeof(unsigned)*(m+1));
  pmat->N_tot = (float*)malloc(sizeof(float)*m);
  nums = i-m;
  pmat->offsets[m] = nums;
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
      i++;
    }
    else{
      //line of iM
      pmat->N_tot[m]  = atof(line);
      pmat->offsets[m++]=i;
    }
    clearline(line);
  }

  fclose(file);

  return nums;
}

void freeMaterialData(struct matdata* pdata){
  free(pdata->offsets);
  free(pdata->N_tot);
  free(pdata->isotopes);
  free(pdata->densities);   
  free(pdata);
}


void clearline(char *line){
  int j;
  for(j=0;j<MAXLEN;j++)
    line[j]=' ';
}
