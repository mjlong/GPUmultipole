#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include <stdio.h>
#include <string.h>

/*#if defined (__QUICKWC)
#if defined(__CFLOAT)
__constant__ float2 table[LENGTH][LENGTH];
#else
__constant__ double2 table[LENGTH][LENGTH];
#endif
//__constant__ CMPTYPE table[LENGTH*LENGTH*2];
#endif
*/

extern void anyvalue(struct multipoledata*,unsigned,struct matdata*, unsigned, unsigned, unsigned, unsigned, unsigned);
int init_data(char* input, char filenames[][FILENAMELEN]);

void printbless();
int main(int argc, char **argv){
  printbless();

  int numIso,totIso;
//read isotopes
  numIso = count_isotopes(argv[5]);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  isotope_read(argv[5],isotopes);
//read materials
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[6]); 
//move on to device settings
  anyvalue(isotopes,numIso,pmat, totIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));

  return 0;
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
  printf("        .\"\" \'<  `.___\\_<|>_/___.\' >\' \"\".      \n");
  printf("       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |    \n");
  printf("       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /    \n");
  printf("   =====`-.____`.___ \\_____/___.-`___.-'===== \n");
  printf("                     `=---='                  \n");
  printf("\n\n\n");
  printf("   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("        佛祖镇楼                  BUG辟易     \n");
}


