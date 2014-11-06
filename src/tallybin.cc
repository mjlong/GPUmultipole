#include "tallybin.h"

void readbins(double** binedges, char* input){
  FILE* fp = NULL;
  char line[MAXLEN];
  fp = fopen(input,"r");
  unsigned i = 0;
  while(fgets(line,MAXLEN,fp))
    i++;
  fclose(fp);
  *binedges = (double*)malloc(sizeof(double)*i);
  fp = NULL;
  fp = fopen(input,"r");
  for(int j=0;j<i;j++){
    clearline(line);
    fgets(line,MAXLEN,fp);
    (*binedges)[j]=atof(line); 
  }
  fclose(fp);
}
