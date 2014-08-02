#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "global.h"
#include <stdio.h>
#include <string.h>

#define FILENAMELEN 20
#define MAXISOTOPES 10
/*#if defined (__QUICKWC)
#if defined(__CFLOAT)
__constant__ float2 table[LENGTH][LENGTH];
#else
__constant__ double2 table[LENGTH][LENGTH];
#endif
//__constant__ CMPTYPE table[LENGTH*LENGTH*2];
#endif
*/
extern void h5read(struct multipoledata & pole, char filename[]);
extern void anyvalue(struct multipoledata*,unsigned, unsigned, unsigned, unsigned, unsigned);

int main(int argc, char **argv){
  int numIso=-1;
  FILE *fp = NULL;
  char line[MAXISOTOPES][FILENAMELEN];
  struct multipoledata *isotopes;//U238data;
  fp = fopen(argv[5],"r");

  if (fp == NULL) {
    fprintf(stderr, "Can't open input file in.list!\n");
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
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  for(i=0;i<numIso;i++)  
    h5read(isotopes[i],line[i]);

  anyvalue(isotopes,numIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  free(isotopes);
  fclose(fp);
  return 0;
}
