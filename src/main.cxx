#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include <stdio.h>
#include <string.h>

#include <optix.h>
#include <cudpp.h>
#include <cudpp_config.h>
//extern void anyvalue(struct multipoledata* data, unsigned numIsos, struct matdata* pmat, unsigned totIsos, unsigned setgridx, unsigned setblockx, unsigned num_src, unsigned devstep);
void printbless();
int main(int argc, char **argv){
  printbless();

//create context
  RTcontext context;
  rtContextCreate(&context);
  int id=0;
  rtContextSetDevices(context, 1, &id);


//CUDPP
//Initialize CUDPP
  CUDPPHandle theCudpp;
  cudppCreate(&theCudpp);
  CUDPPConfiguration config;
  config.datatype = CUDPP_DOUBLE;
  config.algorithm = CUDPP_SORT_RADIX;
  config.options=CUDPP_OPTION_KEY_VALUE_PAIRS;

  CUDPPHandle sortplan = 0;
  CUDPPResult res = cudppPlan(theCudpp, &sortplan, config, atoi(argv[1])*atoi(argv[2]), 1, 0);

  if (CUDPP_SUCCESS != res)
  {
      printf("Error creating CUDPPPlan\n");
      exit(-1);
  }

  res = cudppDestroyPlan(sortplan);
  if (CUDPP_SUCCESS != res)
  {
      printf("Error destroying CUDPPPlan\n");
      exit(-1);
  }
  // shut down the CUDPP library
  cudppDestroy(theCudpp);


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
//  anyvalue( isotopes,numIso,pmat, totIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));


  rtContextDestroy(context); 
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


