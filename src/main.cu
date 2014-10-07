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

#include <cuda.h>
#include "gpuerrchk.h"
extern void h5read(struct multipoledata & pole, char filename[]);
extern void anyvalue(struct multipoledata*,unsigned, unsigned, unsigned, unsigned, unsigned);
int init_data(char* input, char filenames[][FILENAMELEN]);

void printbless();
int main(int argc, char **argv){
  printbless();

  int numIso;
  char filenames[MAXISOTOPES][FILENAMELEN];
  numIso = init_data(argv[5],filenames);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  for(int i=0;i<numIso;i++)  
    h5read(isotopes[i],filenames[i]);

  anyvalue(isotopes,numIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));

  return 0;
}

int init_data(char* input, char line[][FILENAMELEN]){
  int numIso=-1;
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
  return numIso;
}









void printbless(){
  printf("                     _oo0oo_                  \n");  
  printf("                    o8888888o                 \n");
  printf("                    88\" . \"88                 \n");
  printf("                    (| -_- |)                 \n"); 
  printf("                    0\\  =  /0                 \n");
  printf("                  ___/`---'\\___               \n");
  printf("                .' \\\\|     |// '.             \n");
  printf("               / \\\\|||  :  |||// \\            \n");
  printf("              / _||||| -:- |||||- \\           \n");
  printf("             |   | \\\\\\  -  /// |   |          \n");  
  printf("             | \\_|  ''\\---/''  |_/ |          \n");
  printf("             \\  .-\\__  '-'  ___/-. /          \n");
  printf("           ___'. .'  /--.--\\  `. .'___        \n");
  printf("        ."" '<  `.___\\_<|>_/___.' >' "".      \n");
  printf("       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |    \n");
  printf("       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /    \n");
  printf("   =====`-.____`.___ \\_____/___.-`___.-'===== \n");
  printf("                     `=---='                  \n");
  printf("\n\n\n");
  printf("   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("        佛祖镇楼                  BUG辟易     \n");
}


