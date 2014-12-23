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
  unsigned num_batch;
  unsigned num_src;
  num_batch = atoi(argv[1]);
  num_src = atoi(argv[2]);
//============================================================ 
//=============simulation memory allocation===================
//============================================================
  initialize_device();
  MemStruct HostMem;
  unsigned num_bin = readbins(&(HostMem.tallybins),"tallybins")-1;
  initialize_memory(&HostMem, num_bin);
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
  numIso = count_isotopes(argv[3]);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  isotope_read(argv[3],isotopes);
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
  totIso=matread(pmat,argv[4]); 
//copy host material setting to device
  material mat(pmat, totIso);
#if defined(__XS_GPU) 
  unsigned maxiso = 0;
  for(int imat=0;imat<pmat->numMat;imat++){
    int numiso;
    printf("imat=%d:",imat);
    for(int iiso=pmat->offsets[imat];iiso<pmat->offsets[imat+1];iiso++){
       printf(" %d,",pmat->isotopes[iiso]);
    }
    numiso = pmat->offsets[imat+1] - pmat->offsets[imat];
    if(numiso>maxiso)
      maxiso=numiso;
    printf("\n");
  } 
  printf("maxiso=%d\n",maxiso);
  int *iS_h = (int*)malloc(sizeof(int)*maxiso); 
  CMPTYPE *sigTs = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
  CMPTYPE *sigAs = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
  CMPTYPE *sigFs = (CMPTYPE*)malloc(sizeof(CMPTYPE)*maxiso);
#endif
//============================================================ 
//===============main simulation body=========================
//============================================================
clock_t clock_start, clock_end;
float time_elapsed = 0.f;
unsigned active,ibatch,ihistory;
double energy,rnd;
srand(0);
HostMem.num_terminated_neutrons=0;
clock_start = clock();

ibatch = 0;
while(ibatch<num_batch){
ihistory = 0;
while(ihistory<num_src){
active = 1u;
energy = STARTENE;
while(active){
  HostMem.spectrum[search_bin(energy,HostMem.tallybins)]+=1;
  rnd = rand()/(double)RAND_MAX;  
  energy = energy*rnd;
  active = energy>ENDENERG;
  HostMem.tally.cnt+=1;
}
  HostMem.num_terminated_neutrons+=1;
  ihistory+=1;
}
  ibatch+=1;
}
clock_end   = clock();
time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
printf("[time], active cycles costs %f ms\/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons);
print_results(num_src, num_bin, HostMem, time_elapsed);
 
//============================================================ 
//=============simulation shut down===========================
//============================================================
  free(iS_h);
  free(sigTs);
  free(sigAs);
  free(sigFs); 
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


