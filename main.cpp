#include "h5_rdwt.h"
#include "multipole.h"
#define NUM 10
int main(){
  char h5filename[] = "092238.h5";
  double sigT, sigA, sigF;
  double sibT, sibA, sibF;
  double e0=10.0,e1=100000.0,rnd, energy;
  double T = 900.0; 
  multipole U238;
  FILE *file = fopen("results.txt","w");
  

  h5read(U238, h5filename);
  srand(0);
  for(int i=0;i<NUM;i++){
    //    rnd = rand()/(double)RAND_MAX;
    //    energy = e0 + rnd*(e1-e0);
    energy = (i+1.0)*1.63;
    U238.xs_eval_fast(energy,sqrt(KB*T),sigT,sigA,sigF);
    U238.xs_eval_fast(energy           ,sibT,sibA,sibF);
    fprintf(file, "%8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n", 
	    energy, sigT, sibT, sigA, sibA, sigF, sibF);
  }
  fclose(file);

  return 0;
}
