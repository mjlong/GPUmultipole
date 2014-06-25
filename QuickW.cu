#include "QuickW.h"
#if defined(__QUICKWT)
//texture<float2> tex_wtable;
/*static __inline__ __device__ CComplex<float> texfetch_complex8(texture<float2> t, int i){
  float2 v = tex1Dfetch(t,i);
  return CComplex<float>(v.x, v.y);
  }
*/
#if defined(__CFLOAT)
extern texture<float2,2> tex_wtable;
static __inline__ __device__ CComplex<float> texfetch_complex(texture<float2,2> t, int i, int j){
  float2 v = tex2D(t, i, j);
  return CComplex<float>(v.x, v.y);
}
#else
extern texture<int4,2> tex_wtable;
static __inline__ __device__ CComplex<double> texfetch_complex(texture<int4,2> t, int i, int j){
  int4 v = tex2D(t, i, j);
  return CComplex<double>(__hiloint2double(v.y, v.x),__hiloint2double(v.w,v.z));
}
#endif
#endif

#if defined(__QUICKWC)
#if defined(__CFLOAT)
extern __constant__ float2 table[LENGTH*LENGTH];
#else
extern __constant__ double2 table[LENGTH*LENGTH];
#endif
//extern __constant__ CMPTYPE table[LENGTH*LENGTH*2];
#endif

__device__ CMPTYPE b = 0.275255128608410950901357962647054304017026259671664935783653;
__device__ CMPTYPE d = 2.724744871391589049098642037352945695982973740328335064216346;
__device__ CMPTYPE a = 0.512424224754768462984202823134979415014943561548661637413182;
__device__ CMPTYPE c = 0.051765358792987823963876628425793170829107067780337219430904;

/*===============================================================================  
 INITIALIZE_W_TABULATED calculates the Faddeeva function on a 62 x 62 grid       
 using libcerf which is based on Faddeeva package                                                 
 (http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package). The implementation                  
 has accuracy of at least 13 significant digits.                                                 
===============================================================================*/     

/*
__device__ void initialize_w_tabulated(CComplex* w_tabulated){
  int i,j;
  double x,y;
  CComplex z;
  for(i=0;i<LENGTH;i++){
    y = WIDTH*(i-1);
    for(j=0;j<LENGTH;j++){
      x = WIDTH*(j-1);
      z = CComplex(x,y);
      w_tabulated[i*LENGTH+j] = Faddeeva::w(z);
    }
  }
  return;
}
*/

__device__ void fill_w_tabulated(CComplex<CMPTYPE>* w_tabulated, int id){
  double x,y;
  CComplex<double> z;
  y = WIDTH*(id/LENGTH-1);
  x = WIDTH*(id%LENGTH-1);
  z = CComplex<double>(x,y);
#if defined(__CFLOAT)
  z=Faddeeva::w(z);
  w_tabulated[id] = CComplex<float>((float)real(z),(float)imag(z));
#else
  w_tabulated[id] = Faddeeva::w(z);
#endif
  return;
}


/*===============================================================================                   
 W_FUNCTION calculates the Faddeeva function, also known as the complex                
 probability integral, for complex arguments. For |z| < 6, it uses a six-point 
 interpolation scheme based on pre-tabulated data that is accurate to          
 O(10^-3). For |z| > 6, it uses a three-term asymptotic approximation that is                 
 accurate to O(10^-6).                           
===============================================================================*/ 
#if defined(__QUICKWT) || defined(__QUICKWC)
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
  CComplex<CMPTYPE> w;
  if(abs(Norm(z)) < 6.0){
    CMPTYPE p = 10.0*abs(real(z));
    CMPTYPE q = 10.0*imag(z);
    int     l = (int)p + 1;
    int     m = (int)q + 1;
    p = p - (int)p;
    q = q - (int)q;
    CMPTYPE pp = p*p;
    CMPTYPE qq = q*q;
    CMPTYPE pq = p*q;
    w =  
#if defined(__QUICKWT)
      (CMPTYPE)0.5*(qq - q)        *texfetch_complex(tex_wtable, m-1, l  ) + 
      (CMPTYPE)0.5*(pp - p)        *texfetch_complex(tex_wtable, m  , l-1) +
      (CMPTYPE)(1.0 + pq - pp - qq)*texfetch_complex(tex_wtable, m  , l  ) +
      (CMPTYPE)(0.5*(pp + p) - pq) *texfetch_complex(tex_wtable, m  , l+1) +
      (CMPTYPE)(0.5*(qq + q) - pq) *texfetch_complex(tex_wtable, m+1, l  ) +
      (CMPTYPE) pq                 *texfetch_complex(tex_wtable, m+1, l+1);
#else    // __QUICKWC
      (CMPTYPE)0.5*(qq - q)        *CComplex<CMPTYPE>(table[((m-1)*LENGTH+l)*2],table[((m-1)*LENGTH+l)*2+1]),
      (CMPTYPE)0.5*(pp - p)        *CComplex<CMPTYPE>(table[(m*LENGTH + l-1)*2],table[(m*LENGTH + l-1)*2+1]),
      (CMPTYPE)(1.0 + pq - pp - qq)*CComplex<CMPTYPE>(table[(m*LENGTH + l  )*2],table[(m*LENGTH + l  )*2+1]),
      (CMPTYPE)(0.5*(pp + p) - pq) *CComplex<CMPTYPE>(table[(m*LENGTH + l+1)*2],table[(m*LENGTH + l+1)*2+1]),
      (CMPTYPE)(0.5*(qq + q) - pq) *CComplex<CMPTYPE>(table[((m+1)*LENGTH+l)*2],table[((m+1)*LENGTH+l)*2+1]),
      (CMPTYPE) pq                 *CComplex<CMPTYPE>(table[((m+1)*LENGTH+l+1)*2],table[((m+1)*LENGTH+l+1)*2]);
#endif
    if(real(z)<0)
      w = Conjugate(w);
  }
  else
    w = ONEI * z * (a/(z*z - b) + c/(z*z - d));
  return w;
  
}
#else //__QUICKWG
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z, CComplex<CMPTYPE>* w_tabulated){
  CMPTYPE  p;           // interpolation factor on real axis                                   
  CMPTYPE  q;           // interpolation factor on imaginary axis                                  
  CMPTYPE  pp, qq, pq;  // products of p and q                                         

  int l;               //interpolation index for real axis
  int m;               //interpolation index for imaginary axis
  
  CComplex<CMPTYPE> w;
  
  if(abs(Norm(z)) < 6.0){
    // Use interpolation for |z| < 6. The interpolation scheme uses a bivariate         
    // six-point quadrature described in Abramowitz and Stegun 25.2.67. This          
    // interpolation is accurate to O(h^3) = O(10^-3).                           
    //                          
    //     l-1  l  l+1            
    // m+1      +   +                                   
    //          |                         
    // m    +---+---+                         
    //          |                   
    // m-1      +        

    // Determine indices on grid for interpolation and interpolation factors --
    // note that in previous implementations it was necessary to add/subtract
    // two in places because of the indexing on the tabulated function. Because
    // w_tabulated is indexed from -1 to 60, we don't need to do that here

    p = 10.0*abs(real(z));
    q = 10.0*imag(z);
    l = (int)p;
    m = (int)q;
    p = p - l;
    q = q - m;

    //Calculate products
    pp = p*p;
    qq = q*q;
    pq = p*q;

    // Use six-point interpolation to calculate real and imaginary parts
    l++;
    m++;
    w =  
      /*a_b*w_tabulated[(m-1)*LENGTH+l] + 
      a_l*w_tabulated[m*LENGTH + l-1] +
      a_c*w_tabulated[m*LENGTH + l  ] +
      a_r*w_tabulated[m*LENGTH + l+1] +
      a_t*w_tabulated[(m+1)*LENGTH+l] +
      pq *w_tabulated[(m+1)*LENGTH+l+1];*/
      (CMPTYPE)0.5*(qq - q)        *w_tabulated[(m-1)*LENGTH+l] + 
      (CMPTYPE)0.5*(pp - p)        *w_tabulated[m*LENGTH + l-1] +
      (CMPTYPE)(1.0 + pq - pp - qq)*w_tabulated[m*LENGTH + l  ] +
      (CMPTYPE)(0.5*(pp + p) - pq) *w_tabulated[m*LENGTH + l+1] +
      (CMPTYPE)(0.5*(qq + q) - pq) *w_tabulated[(m+1)*LENGTH+l] +
      (CMPTYPE) pq                 *w_tabulated[(m+1)*LENGTH+l+1];

    if(real(z)<0) 
      w = Conjugate(w);
  }
  else
    w = ONEI * z * (a/(z*z - b) + c/(z*z - d));

  return w;
}
#endif

