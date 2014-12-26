#include "CPUComplex.h"
#include "multipole_data.h"
#include "material_data.h"
#include "tallybin.h"
#include <stdio.h>
#include <string.h>
#define KB (8.617342E-5)

//#include <cudpp.h>
//#include <cudpp_config.h>
#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"
#include "simulation.h"

#if defined(__XS_GPU)
#include "multipole.h"
#include "material.h"
#endif

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
#if defined(__XS_GPU)
#if defined(__QUICKWG)
  multipole mp_para(isotopes, numIso, wtable);
#else
  multipole mp_para(isotopes, numIso);
#endif 
#endif
//============================================================ 
//=======Read Materials([isotope, density] pairs)=============
//============================================================
//read from text setting file to host memory 
  struct matdata *pmat=(struct matdata*)malloc(sizeof(struct matdata));
  totIso=matread(pmat,argv[4]); 
//copy host material setting to device
  material mat(pmat, totIso);
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
#if defined(__XS_GPU) 
  unsigned *iS_h, *iS_d; 
  CMPTYPE *sigTs_h,*sigAs_h,*sigFs_h,*sigTs_d,*sigAs_d,*sigFs_d;
  allocate_buffer(maxiso,&iS_d,&sigTs_h,&sigAs_h,&sigFs_h,&sigTs_d,&sigAs_d,&sigFs_d); 
#else
  CMPTYPE sigT, sigA, sigF;
#endif
//============================================================ 
//===============main simulation body=========================
//============================================================
clock_t clock_start, clock_end;
float time_elapsed = 0.f;
unsigned active,ibatch,ihistory;
unsigned imat,iso,numiso;
double energy,sigTsum,sigAsum,sigFsum,rnd;
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
  imat = 0;//fix it to 0 for the moment
#if defined(__XS_GPU)
  iS_h = pmat->isotopes+pmat->offsets[imat];
  numiso = pmat->offsets[imat+1] - pmat->offsets[imat];
  eval_xs(mp_para, iS_h,iS_d, numiso, energy,sqrt(300.0*KB), 
          sigTs_h, sigAs_h, sigFs_h, sigTs_d, sigAs_d, sigFs_d);
#endif
  HostMem.spectrum[search_bin(energy,HostMem.tallybins)]+=1;
  sigTsum=0;
  sigAsum=0;
  sigFsum=0;
  unsigned iiso=0;
  for(int ii=pmat->offsets[imat];ii<pmat->offsets[imat+1];ii++){
#if defined(__XS_GPU)
    sigTsum += sigTs_h[iiso]*pmat->densities[ii];
    sigAsum += sigAs_h[iiso]*pmat->densities[ii];
    sigFsum += sigFs_h[iiso]*pmat->densities[ii];
    iiso++;
#else
#if defined(__FOURIERW)
    host_xs_eval_fast(isotopes[ii], da,db,energy, sqrt(300.0*KB), 
                      sigT, sigA, sigF);
#else //__MIT
    host_xs_eval_fast(isotopes[ii], energy, sqrt(300.0*KB), 
                      sigT, sigA, sigF);
#endif
    sigTsum += sigT*pmat->densities[ii];
    sigAsum += sigA*pmat->densities[ii];
    sigFsum += sigF*pmat->densities[ii];
#endif
  }
#if defined(__PRINTTRACK__)
  printf("[%2d,%3d] %.14e %.14e %.14e %.14e\n", ibatch,ihistory,energy,sigTsum,sigAsum,sigFsum);
#endif
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
#if defined(__XS_GPU)
  release_buffer(iS_d,sigTs_h,sigAs_h,sigFs_h,sigTs_d,sigAs_d,sigFs_d);
#endif
  free(HostMem.tallybins);
//release host isotope data memory
  freeMultipoleData(numIso,isotopes);
//release host material memory
  freeMaterialData(pmat);

  release_memory(HostMem);
#if defined(__XS_GPU)
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


