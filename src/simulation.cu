#include "simulation.h" 

unsigned search_bin(CMPTYPE energy,float* spectrumbins){
  for(int i=0;i<NUM_BINS;i++){
    if( (spectrumbins[i]>=energy)&&(spectrumbins[i+1]<energy) ) 
      return i;
  }
  return 0;
}

__global__ void device_xs_eval(multipole mp_para, unsigned *iS, CMPTYPE E, CMPTYPE sqrtKT, 
                               CMPTYPE* sigTs, CMPTYPE* sigAs, CMPTYPE* sigFs){
  //use one block whose number of threads is the number of isotopes to evaluate xs
  unsigned id = threadIdx.x;
  mp_para.xs_eval_fast(iS[id], E, sqrtKT, sigTs[id],sigAs[id],sigFs[id]);
}



