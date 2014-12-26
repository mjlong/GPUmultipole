#include "fourierw.hh"

#if defined(__QUICKWF)
CMPTYPE bb = 0.275255128608410950901357962647054304017026259671664935783653;
CMPTYPE dd = 2.724744871391589049098642037352945695982973740328335064216346;
CMPTYPE aa = 0.512424224754768462984202823134979415014943561548661637413182;
CMPTYPE cc = 0.051765358792987823963876628425793170829107067780337219430904;
#endif 

void fill_a(CMPTYPE *a, CMPTYPE *b){
  for(int n=1;n<=M;n++){
    a[n] = exp(-n*n*PI*PI/(taom*taom));
    b[n] = (CMPTYPE)n*n*PI*PI;
  }
}

CPUComplex<CMPTYPE> w_function(CPUComplex<CMPTYPE> z, CMPTYPE*a, CMPTYPE* b){
#if defined(__QUICKWF)
  if(Norm(z) > 6.0){
    return ONEI * z * (aa/(z*z - bb) + cc/(z*z - dd));
    //i know imag(z) must > 0, otherwise mangle or conjugate equality must be used
  }
  else{
#endif
  CPUComplex<CMPTYPE> A = taom*z;
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
  CPUComplex<CMPTYPE> B = 
    (CMPTYPE)0.5*(qq - q)        *CPUComplex<CMPTYPE>(constwtable[(m-1)*LENGTH+l].x,constwtable[(m-1)*LENGTH+l].y)+	  
    (CMPTYPE)0.5*(pp - p)        *CPUComplex<CMPTYPE>(constwtable[m*LENGTH + l-1].x,constwtable[m*LENGTH + l-1].y)+	  
    (CMPTYPE)(1.0 + pq - pp - qq)*CPUComplex<CMPTYPE>(constwtable[m*LENGTH + l  ].x,constwtable[m*LENGTH + l  ].y)+	  
    (CMPTYPE)(0.5*(pp + p) - pq) *CPUComplex<CMPTYPE>(constwtable[m*LENGTH + l+1].x,constwtable[m*LENGTH + l+1].y)+	  
    (CMPTYPE)(0.5*(qq + q) - pq) *CPUComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l].x,constwtable[(m+1)*LENGTH+l].y)+	  
    (CMPTYPE) pq                 *CPUComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l+1].x,constwtable[(m+1)*LENGTH+l+1].y);
#else
  CPUComplex<CMPTYPE> B = exp(ONEI*A);      
#endif
  CPUComplex<CMPTYPE> C = A*A;
  CPUComplex<CMPTYPE> w;
  w = CPUComplex<CMPTYPE>((CMPTYPE)0.0,(CMPTYPE)0.0); 
  for(int n=1;n<M;n++){
    //w = w + exp(-n*n*PI*PI/(taom*taom))*( (CMPTYPE)(((n<<31)>>31)|0x00000001)*B - (CMPTYPE)1.0)/(n*n*PI*PI - C);
    w = w + a[n]*( -B - (CMPTYPE)1.0)/(b[n] - C);
    n++;
    w = w + a[n]*(  B - (CMPTYPE)1.0)/(b[n] - C);
  }
//  w = w + a[M] *( -B - (CMPTYPE)1.0)/(b[M]- C);

  w = (w+w)*A;
  w = w + ((CMPTYPE)1.0-B)/A;
  return w*ONEI;
#if defined(__QUICKWF)
  }
#endif
}

