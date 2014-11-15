#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include "tallybin.h"
#include <stdio.h>
#include <string.h>

#include <optix.h>
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
//release host isotope data memory
  freeMultipoleData(numIso,isotopes);
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
//===============Optix Ray Tracing Context====================
//============================================================
float time_elapsed = 0.f;
clock_t clock_start, clock_end;
clock_start = clock();
  RTcontext context;
  RT_CHECK_ERROR(rtContextCreate(&context));
  int id=0;
  rtContextSetDevices(context, 1, &id);
  float geoPara[6] = {0.48f,0.5f,50.f,1.2f,100.f,100.f};
  //float geoPara[6] = {0.00048f,0.0005f,0.050f,0.0012f,0.100f,0.100f};
                      //r1,  r2,  h/2, p,   t,    H/2
#if defined(__QUICKW)
  initialize_context(context, gridsize, atoi(argv[4]),
                     atoi(argv[5]),atoi(argv[6]), 
                     geoPara, DeviceMem, mp_para,mat, wtable);
#else
  initialize_context(context, gridsize, atoi(argv[4]),
                     atoi(argv[5]),atoi(argv[6]), 
                     geoPara, DeviceMem, mp_para,mat);
#endif
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], initializing context costs %f ms\n", time_elapsed);

//============================================================ 
//===============main simulation body=========================
//============================================================
unsigned active;
active = 1u;
set_ray_tracing_program(context,active);
compile_context(context);
initialize_neutrons(gridx, blockx, DeviceMem); 
clock_start = clock();
while(active){
  //since transport_neutrons() surrects all neutrons, rtLaunch always works full load, no need to sort here
  RT_CHECK_ERROR(rtContextLaunch1D(context, 0, gridsize));
  active = count_neutrons(gridx,blockx,DeviceMem,HostMem,num_src);
}
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], active cycles costs %f ms/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);

set_ray_tracing_program(context,active);
compile_context(context);
//while(active){
  //since transport_neutrons() surrects all neutrons, rtLaunch always works full load, no need to sort here
  RT_CHECK_ERROR(rtContextLaunch1D(context, 0, gridsize));
  active = count_neutrons(gridx,blockx,DeviceMem,HostMem,num_src);
//}
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], all cycles costs %f ms/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);

print_results(gridx, blockx, num_src, num_bin, DeviceMem, HostMem, time_elapsed);
 
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
// destroy the optix ray tracing context
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


