#include "h5_rdwt.h"
#include "multipole.h"
#define NUM 10

__global__ void evaluate(multipole pole, double energy, double sqrKT, double *sigT, double *sigA, double *sigF){
  U238.xs_eval_fast(energy,sqrt(KB*T),sigT,sigA,sigF);
}


int main(){
  size_t size = 2*sizeof(double);
  char h5filename[] = "092238.h5";
  double sigT, sigA, sigF;
  double *dev_sigT, *dev_sigA, *dev_sigF;
  cudaMalloc((void**)&dev_sigT, size);
  cudaMalloc((void**)&dev_sigA, size);
  cudaMalloc((void**)&dev_sigF, size);
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
    evaluate<<<1,1>>>(energy, sqrt(KB*T), dev_sigT, dev_sigA, dev_sigF);
    cudaMemcpy(&sigT, dev_sigT, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&sigA, dev_sigA, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&sigF, dev_sigF, size, cudaMemcpyDeviceToHost);
    fprintf(file, "%8.4f %8.5e %8.5e %8.5e \n", 
	           energy, sigT, sigA, sigF);
  }
  fclose(file);
  cudaFree(dev_sigT);
  cudaFree(dev_sigA);
  cudaFree(dev_sigF);
  return 0;
}
