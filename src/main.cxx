#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include "tallybin.h"
#include <stdio.h>
#include <string.h>

//#include <cudpp.h>
//#include <cudpp_config.h>
#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "multipole.h"
#include "material.h"
#include "simulation.h"

#include <time.h>
void printbless();
int main(int argc, char **argv){
//printbless();
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
  MemStruct HostMem;
  unsigned num_bin = readbins(&(HostMem.tallybins),"tallybins")-1;
  initialize_memory(&HostMem, num_bin, gridx,blockx);
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
//============================================================ 
//=======Read Materials([isotope, density] pairs)=============
//============================================================
//read from text setting file to host memory 
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[8]); 
//copy host material setting to device
  material mat(pmat, totIso);
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
clock_start = clock();
//energy = STARTENE;
while(active){
  active=0;  
}
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], active cycles costs %f ms\/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);
print_results(gridx, blockx, num_src, num_bin, HostMem, time_elapsed);
 
//============================================================ 
//=============simulation shut down===========================
//============================================================
  free(HostMem.tallybins);
//release host isotope data memory
  freeMultipoleData(numIso,isotopes);
//release host material memory
  freeMaterialData(pmat);

  release_memory(HostMem);
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
  return 0;
}







#define BLESS " "//"[缪]"
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


