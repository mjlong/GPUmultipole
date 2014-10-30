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
#include "material.h"
#include "optixmain.h"

#include <time.h>
void printbless();
int main(int argc, char **argv){
// 
//calculation dimension
//
  unsigned gridx, blockx, gridsize;
  unsigned num_src, devstep;
  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
  num_src = atoi(argv[3]);
  devstep = atoi(argv[4]);

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
// construct coefficients a[n] for fourier expansion w
#if defined(__FOURIERW)
  CMPTYPE *da;
  CMPTYPE *db;
  fill_wtables(&da,&db);
#endif
// fill exp(z) table for fourierw
#if defined(__INTERPEXP)
  CComplex<CMPTYPE> *exptable;
  fill_wtables(&exptable);
#endif
// fill w function table for Quick W
#if defined(__QUICKW)
  CComplex<CMPTYPE> *wtable;
  fill_wtables(&wtable);
#endif

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



  int numIso,totIso;
//read isotopes
  numIso = count_isotopes(argv[5]);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  isotope_read(argv[5],isotopes);
#if defined(__QUICKWG)
  multipole U238(isotopes, numIso, wtable);
#else
  multipole U238(isotopes, numIso);
#endif 
  freeMultipoleData(numIso,isotopes);

//read materials
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[6]); 
  material mat(pmat, totIso);
  freeMaterialData(pmat);

//move on to device settings
//  anyvalue( isotopes,numIso,pmat, totIso, gridx,blockx,atoi(argv[3]),atoi(argv[4]),cnt,blockcnt, hostarray,devicearray, HostMem,DeviceMem);

//
//main simulation body
//
clock_t clock_start, clock_end;
float time_elapsed = 0.f;
unsigned active;
#if defined(__PROCESS) //|| defined(__TRACK)
  active = 0u;
#else
  active = 1u;
#endif
initialize_neutrons(gridx, blockx, DeviceMem); 
clock_start = clock();

while(active){
  start_neutrons(gridx, blockx, numIso, U238, devicearray, DeviceMem, num_src, devstep);
  active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
}
  remain_neutrons(gridx, blockx,numIso, U238, devicearray, DeviceMem);
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
print_results(gridx, blockx, num_src, devstep, DeviceMem, HostMem, hostarray, devicearray, blockcnt,cnt, time_elapsed);
 
  release_memory(DeviceMem, HostMem, devicearray, hostarray, cnt, blockcnt);
  U238.release_pointer();
  mat.release_pointer();

#if defined(__FOURIERW)
  release_wtables(da,db);
#endif
#if defined(__INTERPEXP)
  release_wtables(exptable);
#endif
#if defined(__QUICKW)
  release_wtables(wtable);
#endif
  res = cudppDestroyPlan(sortplan);
  if (CUDPP_SUCCESS != res)
  {
      printf("Error destroying CUDPPPlan\n");
      exit(-1);
  }
  // shut down the CUDPP library
  cudppDestroy(theCudpp);

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


