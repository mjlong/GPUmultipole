#include "simulation.h" 

unsigned search_bin(CMPTYPE energy,float* spectrumbins){
  for(int i=0;i<NUM_BINS;i++){
    if( (spectrumbins[i]>=energy)&&(spectrumbins[i+1]<energy) ) 
      return i;
  }
  return 0;
}
