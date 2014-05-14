#include "h5_rdwt.h"
#include "multipole.h"
#define NUM 10
int main(){
  char h5filename[] = "092238.h5";
  double sigT, sigA, sigF;
  double sibT, sibA, sibF;
  double e0=1.0,e1=2000.0,rnd, energy;
  double T = 300.0; 
  unsigned *counts;
  multipole U238;
  FILE *file = fopen("results.txt","w");
  
  int i,j;
  counts = (unsigned*)malloc(sizeof(unsigned)*2*N*2*N);
  for(i=0;i<2*N*2*N;i++)
      counts[i]=0;
  h5read(U238, h5filename);
  srand(0);
  e1 = 1.95093e4;
  energy = e1;
  e0 = e1;//note: this extra line disables all loops
  while(energy>=e0){
    rnd = rand()/(double)RAND_MAX;
    //energy = 20000.0; //(i+1.0)*1.63;//20000.0;//
    U238.xs_eval_fast(energy,sqrt(KB*T),sigT,sigA,sigF,counts,2*N);
    printf("xs:%g,%g,%g\n",sigT,sigA,sigF);
    energy = energy*rnd;
  }

  /*  for(i=0;i<2*N;i++){
    for(j=0;j<2*N;j++)
      fprintf(file,"%4d",counts[i*2*N+j]);
    fprintf(file,"\n");
    }*/
  free(counts);
  fclose(file);
  return 0;
}
/*
    U238.xs_eval_fast(energy           ,sibT,sibA,sibF);
    fprintf(file, "%8.4f %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n", 
	    energy, sigT, sibT, sigA, sibA, sigF, sibF);
*/
