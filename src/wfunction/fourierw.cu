#include "fourierw.h"

extern __constant__ CMPTYPE a[];
extern __constant__ CMPTYPE b[];

#if defined(__QUICKWF)
__constant__ CMPTYPE bb = 0.275255128608410950901357962647054304017026259671664935783653;
__constant__ CMPTYPE dd = 2.724744871391589049098642037352945695982973740328335064216346;
__constant__ CMPTYPE aa = 0.512424224754768462984202823134979415014943561548661637413182;
__constant__ CMPTYPE cc = 0.051765358792987823963876628425793170829107067780337219430904;
#endif 

#if defined(__INTERPEXP)
extern __constant__ CMPTYPE2 constwtable[];
__global__ void fill_exp_table(CComplex<CMPTYPE>* exptable){
  int id = blockDim.x*blockIdx.x + threadIdx.x;
  CMPTYPE x,y;
  CComplex<CMPTYPE> z;
  y = WIDTH*(id/LENGTH-1);
  x = WIDTH*(id%LENGTH-1);
  z = CComplex<CMPTYPE>(x,y);
  exptable[id] = exp(z*taom*ONEI);
  return;
}
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
  CComplex<CMPTYPE> A = taom*z;
#if defined(__INTERPEXP)
  CMPTYPE p = 10.0*abs(real(z));
  CMPTYPE q = 10.0*imag(z);
  int     l = (int)p + 1;
  int     m = (int)q + 1;
  p = p - (int)p;
  q = q - (int)q;
  CMPTYPE pp = p*p;
  CMPTYPE qq = q*q;
  CMPTYPE pq = p*q;
  CComplex<CMPTYPE> B = 
    (CMPTYPE)0.5*(qq - q)        *CComplex<CMPTYPE>(constwtable[(m-1)*LENGTH+l].x,constwtable[(m-1)*LENGTH+l].y)+	  
    (CMPTYPE)0.5*(pp - p)        *CComplex<CMPTYPE>(constwtable[m*LENGTH + l-1].x,constwtable[m*LENGTH + l-1].y)+	  
    (CMPTYPE)(1.0 + pq - pp - qq)*CComplex<CMPTYPE>(constwtable[m*LENGTH + l  ].x,constwtable[m*LENGTH + l  ].y)+	  
    (CMPTYPE)(0.5*(pp + p) - pq) *CComplex<CMPTYPE>(constwtable[m*LENGTH + l+1].x,constwtable[m*LENGTH + l+1].y)+	  
    (CMPTYPE)(0.5*(qq + q) - pq) *CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l].x,constwtable[(m+1)*LENGTH+l].y)+	  
    (CMPTYPE) pq                 *CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l+1].x,constwtable[(m+1)*LENGTH+l+1].y);
#else
  CComplex<CMPTYPE> B = exp(ONEI*A);      
#endif
  CComplex<CMPTYPE> C = A*A;
  CComplex<CMPTYPE> w;
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

