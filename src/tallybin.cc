#include "tallybin.h"

unsigned readbins(float** binedges, char* input){
  FILE* fp = NULL;
  char line[MAXLEN];
  fp = fopen(input,"r");
  unsigned i = 0;
  while(fgets(line,MAXLEN,fp))
    i++;
  fclose(fp);
  *binedges = (float*)malloc(sizeof(float)*i);
  fp = NULL;
  fp = fopen(input,"r");
  for(int j=0;j<i;j++){
    clearline(line);
    fgets(line,MAXLEN,fp);
    (*binedges)[j]=atof(line); 
  }
  fclose(fp);
  return i;
}
