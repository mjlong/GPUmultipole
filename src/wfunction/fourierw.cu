#include "fourierw.h"

extern __constant__ CMPTYPE a[];
extern __constant__ CMPTYPE b[];

#if defined(__QUICKWF)
__constant__ CMPTYPE bb = 0.275255128608410950901357962647054304017026259671664935783653;
__constant__ CMPTYPE dd = 2.724744871391589049098642037352945695982973740328335064216346;
__constant__ CMPTYPE aa = 0.512424224754768462984202823134979415014943561548661637413182;
__constant__ CMPTYPE cc = 0.051765358792987823963876628425793170829107067780337219430904;
#endif 

__global__ void fill_a(CMPTYPE *a, CMPTYPE *b){
  // since n is supposed to be 1,2,3,...,23
  // by default dimGrid=(1,1,1)
  int n = threadIdx.x+1;
  a[n] = exp(-n*n*PI*PI/(taom*taom));
  b[n] = (CMPTYPE)n*n*PI*PI;
}

__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
#if defined(__QUICKWF)
  if(Norm(z) > 6.0){
    return ONEI * z * (aa/(z*z - bb) + cc/(z*z - dd));
    //i know imag(z) must > 0, otherwise mangle or conjugate equality must be used
  }
  else{
#endif
  CComplex<CMPTYPE> w;
  CComplex<CMPTYPE> A = taom*z;
  CComplex<CMPTYPE> B = exp(ONEI*A);      
  CComplex<CMPTYPE> C = A*A;
  w = CComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0); 
  for(int n=1;n<M;n++){
    //w = w + exp(-n*n*PI*PI/(taom*taom))*( (CMPTYPE)(((n<<31)>>31)|0x00000001)*B - (CMPTYPE)1.0)/(n*n*PI*PI - C);
    w = w + a[n]*( -B - (CMPTYPE)1.0)/(b[n] - C);
    n++;
    w = w + a[n]*(  B - (CMPTYPE)1.0)/(b[n] - C);
  }
  w = w + a[M] *( -B - (CMPTYPE)1.0)/(b[M]- C);

  w = (w+w)*A;
  w = w + ((CMPTYPE)1.0-B)/A;
  return w*ONEI;
#if defined(__QUICKWF)
  }
#endif
}

