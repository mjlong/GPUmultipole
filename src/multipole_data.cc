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
