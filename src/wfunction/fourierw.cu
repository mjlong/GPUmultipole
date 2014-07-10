#include "fourierw.h"

extern __constant__ CMPTYPE a[];

__global__ void fill_a(CMPTYPE *a){
  // since n is supposed to be 1,2,3,...,23
  // by default dimGrid=(1,1,1)
  int n = threadIdx.x+1;
  a[n] = exp(-n*n*PI*PI/(taom*taom));
}

__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
  CComplex<CMPTYPE> w;
  CComplex<CMPTYPE> A = taom*z;
  CComplex<CMPTYPE> B = exp(ONEI*A);
  CComplex<CMPTYPE> C = A*A;
  w = CComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0); 
  for(int n=1;n<M;n++){
    //w = w + exp(-n*n*PI*PI/(taom*taom))*( (CMPTYPE)(((n<<31)>>31)|0x00000001)*B - (CMPTYPE)1.0)/(n*n*PI*PI - C);
    w = w + a[n]*( -B - (CMPTYPE)1.0)/((CMPTYPE)(n*n*PI*PI) - C);
    n++;
    w = w + a[n]*(  B - (CMPTYPE)1.0)/((CMPTYPE)(n*n*PI*PI) - C);
  }
  w = w + a[23] *( -B - (CMPTYPE)1.0)/((CMPTYPE)(23*23*PI*PI) - C);

  w = w*2.0*A;
  w = w + ((CMPTYPE)1.0-B)/(taom*z);
  w = ONEI*w;
  return w;
}

