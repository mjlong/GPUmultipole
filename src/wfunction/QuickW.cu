#include "QuickW.h"
#if defined(__QUICKWT)
#if defined(__CFLOAT)
texture<float2> tex_wtable;
static __inline__ __device__ CComplex<float> texfetch_complex(texture<float2> t, int i){
  float2 v = tex1Dfetch(t,i);
  return CComplex<float>(v.x, v.y);
  }
#else
texture<int4> tex_wtable;
static __inline__ __device__ CComplex<double> texfetch_complex(texture<int4> t, int i){
  int4 v = tex1Dfetch(t,i);
  return CComplex<double>(__hiloint2double(v.y, v.x),__hiloint2double(v.w,v.z));
}
#endif
/*
#if defined(__CFLOAT)
texture<float2,2> tex_wtable;
static __inline__ __device__ CComplex<float> texfetch_complex(texture<float2,2> t, int i, int j){
  float2 v = tex2D(t, i, j);
  return CComplex<float>(v.x, v.y);
}
#else
texture<int4,2> tex_wtable;
static __inline__ __device__ CComplex<double> texfetch_complex(texture<int4,2> t, int i, int j){
  int4 v = tex2D(t, i, j);
  return CComplex<double>(__hiloint2double(v.y, v.x),__hiloint2double(v.w,v.z));
}
#endif
*/
#endif

#if defined (__QUICKWC)
extern __constant__ CMPTYPE2 constwtable[];
#endif

//TODO: load constants on shared memory
__constant__ CMPTYPE b = 0.275255128608410950901357962647054304017026259671664935783653;
__constant__ CMPTYPE d = 2.724744871391589049098642037352945695982973740328335064216346;
__constant__ CMPTYPE a = 0.512424224754768462984202823134979415014943561548661637413182;
__constant__ CMPTYPE c = 0.051765358792987823963876628425793170829107067780337219430904;

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


__global__ void fill_w_tabulated(CComplex<CMPTYPE>* w_tabulated){
  int id = blockDim.x*blockIdx.x + threadIdx.x;
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

#if defined(__QUICKWT)
__host__ void bindwtable(CComplex<CMPTYPE>* wtable){
  cudaBindTexture(NULL, tex_wtable, wtable, LENGTH*LENGTH*sizeof(CMPTYPE)*2);
  //cudaChannelFormatDesc desc = cudaCreateChannelDesc<CMPTYPE2>();
  //cudaBindTexture2D(NULL, tex_wtable, wtable, desc, LENGTH, LENGTH, sizeof(CMPTYPE)*2*LENGTH);
}

__host__ void unbindwtable(){
  cudaUnbindTexture(tex_wtable);
}
#endif 

/*===============================================================================                   
 W_FUNCTION calculates the Faddeeva function, also known as the complex                
 probability integral, for complex arguments. For |z| < 6, it uses a six-point 
 interpolation scheme based on pre-tabulated data that is accurate to          
 O(10^-3). For |z| > 6, it uses a three-term asymptotic approximation that is                 
 accurate to O(10^-6).                           
===============================================================================*/ 
#if defined(__QUICKWT) 
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
#endif
#if defined(__QUICKWG) 
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z, CComplex<CMPTYPE>* w_tabulated){
#endif 
#if defined(__QUICKWC)
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> z){
#endif
  CComplex<CMPTYPE> w;
  bool lower = (imag(z)<0);
  if(lower) z = Conjugate(z);
  //TODO: if(!lower){}else{...} should be faster but isn't
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
    //**********************************************************************
    //CComplex<CMPTYPE> w1,w2,w3,w4,w5,w6;
    //**********************************************************************
#if defined(__QUICKWT)
    //failed test of binding 1d array to 2d texture
    /*
    w1 = texfetch_complex(tex_wtable, m-1, l  );
    w2 = texfetch_complex(tex_wtable, m  , l-1);
    w3 = texfetch_complex(tex_wtable, m  , l  );
    w4 = texfetch_complex(tex_wtable, m  , l+1);
    w5 = texfetch_complex(tex_wtable, m+1, l  );
    w6 = texfetch_complex(tex_wtable, m+1, l+1);
    */
    /*
    w1 = texfetch_complex(tex_wtable, (m-1)*LENGTH+l);
    w2 = texfetch_complex(tex_wtable, m*LENGTH + l-1);
    w3 = texfetch_complex(tex_wtable, m*LENGTH + l  );
    w4 = texfetch_complex(tex_wtable, m*LENGTH + l+1);
    w5 = texfetch_complex(tex_wtable, (m+1)*LENGTH+l);
    w6 = texfetch_complex(tex_wtable, (m+1)*LENGTH+l+1);
    */
    w =  
      (CMPTYPE)0.5*(qq - q)        *texfetch_complex(tex_wtable,(m-1)*LENGTH+l ) + 
      (CMPTYPE)0.5*(pp - p)        *texfetch_complex(tex_wtable,m*LENGTH + l-1 ) +
      (CMPTYPE)(1.0 + pq - pp - qq)*texfetch_complex(tex_wtable,m*LENGTH + l   ) +
      (CMPTYPE)(0.5*(pp + p) - pq) *texfetch_complex(tex_wtable,m*LENGTH + l+1 ) +
      (CMPTYPE)(0.5*(qq + q) - pq) *texfetch_complex(tex_wtable,(m+1)*LENGTH+l ) +
      (CMPTYPE) pq                 *texfetch_complex(tex_wtable,(m+1)*LENGTH+l+1);
      /*
      (CMPTYPE)0.5*(qq - q)        *texfetch_complex(tex_wtable, m-1, l  ) + 
      (CMPTYPE)0.5*(pp - p)        *texfetch_complex(tex_wtable, m  , l-1) +
      (CMPTYPE)(1.0 + pq - pp - qq)*texfetch_complex(tex_wtable, m  , l  ) +
      (CMPTYPE)(0.5*(pp + p) - pq) *texfetch_complex(tex_wtable, m  , l+1) +
      (CMPTYPE)(0.5*(qq + q) - pq) *texfetch_complex(tex_wtable, m+1, l  ) +
      (CMPTYPE) pq                 *texfetch_complex(tex_wtable, m+1, l+1);
      */
#endif
///////////////////////////////////
#if defined(__QUICKWG)
    /*
    w1 = w_tabulated[(m-1)*LENGTH+l];
    w2 = w_tabulated[m*LENGTH + l-1];
    w3 = w_tabulated[m*LENGTH + l  ];
    w4 = w_tabulated[m*LENGTH + l+1];
    w5 = w_tabulated[(m+1)*LENGTH+l];
    w6 = w_tabulated[(m+1)*LENGTH+l+1];
    */
    w = 
      (CMPTYPE)0.5*(qq - q)        *w_tabulated[(m-1)*LENGTH+l] + 
      (CMPTYPE)0.5*(pp - p)        *w_tabulated[m*LENGTH + l-1] +
      (CMPTYPE)(1.0 + pq - pp - qq)*w_tabulated[m*LENGTH + l  ] +
      (CMPTYPE)(0.5*(pp + p) - pq) *w_tabulated[m*LENGTH + l+1] +
      (CMPTYPE)(0.5*(qq + q) - pq) *w_tabulated[(m+1)*LENGTH+l] +
      (CMPTYPE) pq                 *w_tabulated[(m+1)*LENGTH+l+1];
#endif
////////////////////////////////////////////
#if defined(__QUICKWC)
    /*
    w1 = CComplex<CMPTYPE>(constwtable[(m-1)*LENGTH+l].x,constwtable[(m-1)*LENGTH+l].y);
    w2 = CComplex<CMPTYPE>(constwtable[m*LENGTH + l-1].x,constwtable[m*LENGTH + l-1].y);
    w3 = CComplex<CMPTYPE>(constwtable[m*LENGTH + l  ].x,constwtable[m*LENGTH + l  ].y);
    w4 = CComplex<CMPTYPE>(constwtable[m*LENGTH + l+1].x,constwtable[m*LENGTH + l+1].y);
    w5 = CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l].x,constwtable[(m+1)*LENGTH+l].y);
    w6 = CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l+1].x,constwtable[(m+1)*LENGTH+l+1].y);
    */
    w = 
      (CMPTYPE)0.5*(qq - q)        *CComplex<CMPTYPE>(constwtable[(m-1)*LENGTH+l].x,constwtable[(m-1)*LENGTH+l].y)+	  
      (CMPTYPE)0.5*(pp - p)        *CComplex<CMPTYPE>(constwtable[m*LENGTH + l-1].x,constwtable[m*LENGTH + l-1].y)+	  
      (CMPTYPE)(1.0 + pq - pp - qq)*CComplex<CMPTYPE>(constwtable[m*LENGTH + l  ].x,constwtable[m*LENGTH + l  ].y)+	  
      (CMPTYPE)(0.5*(pp + p) - pq) *CComplex<CMPTYPE>(constwtable[m*LENGTH + l+1].x,constwtable[m*LENGTH + l+1].y)+	  
      (CMPTYPE)(0.5*(qq + q) - pq) *CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l].x,constwtable[(m+1)*LENGTH+l].y)+	  
      (CMPTYPE) pq                 *CComplex<CMPTYPE>(constwtable[(m+1)*LENGTH+l+1].x,constwtable[(m+1)*LENGTH+l+1].y);
#endif
    /*
     w=
      (CMPTYPE)0.5*(qq - q)        *w1+	  
      (CMPTYPE)0.5*(pp - p)        *w2+	  
      (CMPTYPE)(1.0 + pq - pp - qq)*w3+	  
      (CMPTYPE)(0.5*(pp + p) - pq) *w4+	  
      (CMPTYPE)(0.5*(qq + q) - pq) *w5+	  
      (CMPTYPE) pq                 *w6;
    if(blockIdx.x==0 && threadIdx.x==18){
      printf("w1=%20.16e + i*%20.16e\n",real(w1),imag(w1));
      printf("w2=%20.16e + i*%20.16e\n",real(w2),imag(w2));
      printf("w3=%20.16e + i*%20.16e\n",real(w3),imag(w3));
      printf("w4=%20.16e + i*%20.16e\n",real(w4),imag(w4));
      printf("w5=%20.16e + i*%20.16e\n",real(w5),imag(w5));
      printf("w6=%20.16e + i*%20.16e\n",real(w6),imag(w6));
      printf("m=%2d, l=%2d\n", m,l);
      printf("p=%16.12e, q=%16.12e, pp=%16.12e, qq=%16.12e, pq=%16.12e\n",p,q,pp,qq,pq);
    }
    */
    if(real(z)<0) 
      w = Conjugate(w);
  }
  else
    w = ONEI * z * (a/(z*z - b) + c/(z*z - d));
  if(lower){
    //w = Conjugate((CMPTYPE)2.0*exp((CMPTYPE)0.0-z*z)-w);
    w = -Conjugate(w);
#if defined(__PLOT)
    printf("hit\n");
#endif
  }
  return w;
}

