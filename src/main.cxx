#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include "tallybin.h"
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
  printbless();
//============================================================ 
//====================calculation dimension===================
//============================================================
  unsigned gridx, blockx, gridsize;
  unsigned num_src;
  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
  num_src = atoi(argv[3]);
//============================================================ 
//=============simulation memory allocation===================
//============================================================
  initialize_device();
  MemStruct HostMem, DeviceMem;
  unsigned num_bin = readbins(&(HostMem.tallybins),"tallybins")-1;
  initialize_memory(&DeviceMem, &HostMem, num_bin, gridx,blockx);
  free(HostMem.tallybins);
//============================================================ 
//===============Faddeeva tables==============================
//============================================================
//
//===construct coefficients a[n] for fourier expansion w======
//
#if defined(__FOURIERW)
  CMPTYPE *da;
  CMPTYPE *db;
  fill_wtables(&da,&db);
#endif
//
//=============fill exp(z) table for fourierw=================
//
#if defined(__INTERPEXP)
  CComplex<CMPTYPE> *exptable;
  fill_wtables(&exptable);
#endif
//
//==========fill w function table for Quick W=================
//
#if defined(__QUICKW)
  CComplex<CMPTYPE> *wtable;
  fill_wtables(&wtable);
#endif

//============================================================ 
//===============Optix Ray Tracing Context====================
//============================================================
  RTcontext context;
  RT_CHECK_ERROR(rtContextCreate(&context));
  int id=0;
  rtContextSetDevices(context, 1, &id);
  float geoPara[6] = {0.48f,0.5f,50.f,1.2f,100.f,100.f};
  //float geoPara[6] = {0.00048f,0.0005f,0.050f,0.0012f,0.100f,0.100f};
                      //r1,  r2,  h/2, p,   t,    H/2
  initialize_context(context, gridsize, 
                     atoi(argv[5]),atoi(argv[6]), 
                     geoPara, DeviceMem.nInfo);
//============================================================ 
//=============CUDPP Initialization===========================
//============================================================
  CUDPPHandle theCudpp;
  cudppCreate(&theCudpp);
  CUDPPConfiguration config;
  config.datatype = CUDPP_DOUBLE;
  config.algorithm = CUDPP_SORT_RADIX;
  config.options=CUDPP_OPTION_KEY_VALUE_PAIRS;
  config.options=CUDPP_OPTION_BACKWARD;
  CUDPPHandle sortplan = 0;
  CUDPPResult res = cudppPlan(theCudpp, &sortplan, config, gridsize, 1, 0);
  if (CUDPP_SUCCESS != res)
  {
      printf("Error creating CUDPPPlan\n");
      exit(-1);
  }
//============================================================ 
//=============Read Isotopes(multipole data)==================
//============================================================
  int numIso,totIso;
//read from hdf5 file to host memory
  numIso = count_isotopes(argv[7]);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  isotope_read(argv[7],isotopes);
//copy host isotope data to device
#if defined(__QUICKWG)
  multipole mp_para(isotopes, numIso, wtable);
#else
  multipole mp_para(isotopes, numIso);
#endif 
//============================================================ 
//=======Read Materials([isotope, density] pairs)=============
//============================================================
//read from text setting file to host memory 
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[8]); 
//copy host material setting to device
  material mat(pmat, totIso);
//release host material memory
  freeMaterialData(pmat);
//============================================================ 
//===============main simulation body=========================
//============================================================
clock_t clock_start, clock_end;
float time_elapsed = 0.f;
unsigned active;
#if defined(__PROCESS) //|| defined(__TRACK)
  active = 0u;
#else
  active = 1u;
#endif
initialize_neutrons(gridx, blockx, DeviceMem); 
//clock_start = clock();
int after=atoi(argv[4]);
int i=0;
//while(active){
for(i=0;i<after;i++){
  //since transport_neutrons() surrects all neutrons, rtLaunch always works full load, no need to sort here
  RT_CHECK_ERROR(rtContextLaunch1D(context, 0, gridsize));
  //sort key = live*(isotopeID*MAXENERGY+energy)
  sort_prepare(gridx, blockx, DeviceMem, mat);
  cudppRadixSort(sortplan, DeviceMem.nInfo.isoenergy, DeviceMem.nInfo.id, gridsize);
  //                          keys,                   values,             numElements
  //neutrons found leaked in *locate* will not be evaluated 
  start_neutrons(gridx, blockx, mat, mp_para, DeviceMem, num_src,1);
  //besides moving, neutrons terminated is initiated as new 
  active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
  transport_neutrons(gridx, blockx, DeviceMem, mat, active); 
  //if active=1; transport<<<>>> will renew neutrons with live=0
  //if active=0; transport<<<>>> will leave terminated neutrons
  //set active always 1 to make sure number of neutrons simulated exactly equal to num_src
}
//clock_end   = clock();
//time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
//printf("[time], active cycles costs %f ms\/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);
#if defined(__PRINTTRACK__)
unsigned left = count_lives(gridx,blockx,DeviceMem,HostMem);
#else
HostMem.num_terminated_neutrons[0]+=count_lives(gridx,blockx,DeviceMem,HostMem);
#endif
active = 1;

//instead of moving to track the remaining neutrons, i take the array of E's, convert them to z's and play with the array of z's
copyE(HostMem, DeviceMem,gridsize);
for(int i=0;i<gridsize;i++){
  printf("E[%2d]=%.5f\n",i,HostMem.nInfo.energy[i]);
}
CPUComplex<CMPTYPE> *pz;
unsigned numz=generateZ(isotopes[0], sqrt(KB*300.0),HostMem.nInfo.energy,gridsize, &pz);
printf("From %d energies, I have %d complex numbers for Faddeeva evaluation:\n",gridsize,numz);
for(int i=0;i<numz;i++){
  pz[i].output();
}

free(pz);
//release host isotope data memory
  freeMultipoleData(numIso,isotopes);

/*
while(0!=active){
  //about twice sort in one loop
  //1. add extra sort here
  //2. only sort before xs evaluation, allows thread divergence in ray tracing
  sort_prepare(gridx, blockx, DeviceMem, mat);
  cudppRadixSort(sortplan, DeviceMem.nInfo.isoenergy, DeviceMem.nInfo.id, gridsize);
  RT_CHECK_ERROR(rtContextLaunch1D(context, 0, gridsize));

  sort_prepare(gridx, blockx, DeviceMem, mat);
  cudppRadixSort(sortplan, DeviceMem.nInfo.isoenergy, DeviceMem.nInfo.id, gridsize);
  start_neutrons(gridx, blockx, mat, mp_para, DeviceMem, num_src,0);

  active = count_lives(gridx, blockx, DeviceMem, HostMem);
#if defined(__PRINTTRACK__)
  HostMem.num_terminated_neutrons[0]+=left-active;
  left = active;
  printf("[remaining]%d terminated, %d left\n",HostMem.num_terminated_neutrons[0],left);
#endif
  transport_neutrons(gridx, blockx, DeviceMem, mat, 0); 
}
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], active + remain cycles costs %f ms\/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);
print_results(gridx, blockx, num_src, num_bin, DeviceMem, HostMem, time_elapsed);
*/ 
//============================================================ 
//=============simulation shut down===========================
//============================================================
  release_memory(DeviceMem, HostMem);
  mp_para.release_pointer();
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
// destroy the optix ray tracing context
  rtContextDestroy(context); 
  return 0;
}







#define BLESS "[缪]"
void printbless(){
  printf(BLESS"                     _oo0oo_                  \n");  
  printf(BLESS"                    o8888888o                 \n");
  printf(BLESS"                    88\" . \"88                 \n");
  printf(BLESS"                    (| -_- |)                 \n"); 
  printf(BLESS"                    0\\  =  /0                 \n");
  printf(BLESS"                  ___/`---'\\___               \n");
  printf(BLESS"                .' \\\\|     |// '.             \n");
  printf(BLESS"               / \\\\|||  :  |||// \\            \n");
  printf(BLESS"              / _||||| -:- |||||- \\           \n");
  printf(BLESS"             |   | \\\\\\  -  /// |   |          \n");  
  printf(BLESS"             | \\_|  ''\\---/''  |_/ |          \n");
  printf(BLESS"             \\  .-\\__  '-'  ___/-. /          \n");
  printf(BLESS"           ___'. .'  /--.--\\  `. .'___        \n");
  printf(BLESS"        .\"\" \'<  `.___\\_<|>_/___.\' >\' \"\".      \n");
  printf(BLESS"       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |    \n");
  printf(BLESS"       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /    \n");
  printf(BLESS"   =====`-.____`.___ \\_____/___.-`___.-'===== \n");
  printf(BLESS"                     `=---='                  \n");
  printf(BLESS"\n\n\n");
  printf(BLESS"   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf(BLESS"        佛祖镇楼                  BUG辟易     \n");
}


