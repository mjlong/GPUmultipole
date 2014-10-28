#include "material_data.h"

void matread(struct matdata mat, char* filename){
  int i,j;
  int m=0;
  unsigned iS;
  unsigned *isotopes;
  float    *densities;
  float rho;
  FILE *file;
  char line[MAXLEN];
  for(i=0;i<MAXLEN;i++)
    line[i]=' ';
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
  }
  //mat wise data
  mat.offsets = (unsigned*)malloc(sizeof(unsigned)*m);
  mat.N_tot = (float*)malloc(sizeof(float)*m);
  //isotope wise data
  mat.isotopes = (unsigned*)malloc(sizeof(unsigned)*(i-m));
  mat.densities= (float*)malloc(sizeof(float)*(i-m));
  i=0;
  while(fgets(line,MAXLEN,file)){
    for(j=0;j<MAXLEN-1;j++){
      if((' '==line[j]) && (' '!=line[j+1]))
        break;
    }
    printf("line %d: ",++i);
    if(j<MAXLEN-1){
      mat.isotopes[i]  = atoi(line);
      mat.densities[i] = atof(line+j);
    }
    else
      i--;
    for(j=0;j<MAXLEN;j++)
      line[j]=' ';
  }

  fclose(file);

  return;
}


