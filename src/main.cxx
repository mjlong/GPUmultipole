#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include <stdio.h>
#include <string.h>

#include <optix.h>
#include <cudpp.h>
#include <cudpp_config.h>
#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "multipole.h"
#include "optixmain.h"

void printbless();
int main(int argc, char **argv){
// 
//calculation dimension
//
  unsigned gridx, blockx, gridsize;
  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
//
//simulation memory allocation
//
  initialize_device();
  unsigned *cnt, *blockcnt;
  CMPTYPE *hostarray, *devicearray;
  MemStruct HostMem, DeviceMem;
  initialize_memory(&DeviceMem, &HostMem, &devicearray, &hostarray, &cnt, &blockcnt, gridx,blockx);
  printbless();

//
//Faddeeva tables
//

//create context
  RTcontext context;
  rtContextCreate(&context);
  int id=0;
  rtContextSetDevices(context, 1, &id);
  float geoPara[6] = {0.48f,0.5f,50.f,1.2f,100.f,100.f};
  //float geoPara[6] = {0.00048f,0.0005f,0.050f,0.0012f,0.100f,0.100f};
                      //r1,  r2,  h/2, p,   t,    H/2
  tracemain2(context, gridsize, 2, 2, geoPara, DeviceMem.nInfo);
  //rtContextDestroy( context );

//CUDPP
//Initialize CUDPP
  CUDPPHandle theCudpp;
  cudppCreate(&theCudpp);
  CUDPPConfiguration config;
  config.datatype = CUDPP_DOUBLE;
  config.algorithm = CUDPP_SORT_RADIX;
  config.options=CUDPP_OPTION_KEY_VALUE_PAIRS;

  CUDPPHandle sortplan = 0;
  CUDPPResult res = cudppPlan(theCudpp, &sortplan, config, gridsize, 1, 0);

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
  multipole U92238(isotopes, numIso);
  U92238.release_pointer();
//read materials
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[6]); 
//move on to device settings
  anyvalue( isotopes,numIso,pmat, totIso, gridx,blockx,atoi(argv[3]),atoi(argv[4]),cnt,blockcnt, hostarray,devicearray, HostMem,DeviceMem);

 cudppRadixSort(sortplan, DeviceMem.nInfo.isoenergy, DeviceMem.nInfo.id, gridsize);
  printf("tested cudppSort\n");
  release_memory(DeviceMem, HostMem, devicearray, hostarray, cnt, blockcnt);

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


