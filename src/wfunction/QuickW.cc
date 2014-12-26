#include "QuickW.hh"
//TODO: load constants on shared memory
CMPTYPE b = 0.275255128608410950901357962647054304017026259671664935783653;
CMPTYPE d = 2.724744871391589049098642037352945695982973740328335064216346;
CMPTYPE a = 0.512424224754768462984202823134979415014943561548661637413182;
CMPTYPE c = 0.051765358792987823963876628425793170829107067780337219430904;

/*===============================================================================  
 INITIALIZE_W_TABULATED calculates the Faddeeva function on a 62 x 62 grid       
 using libcerf which is based on Faddeeva package                                                 
 (http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package). The implementation                  
 has accuracy of at least 13 significant digits.                                                 
===============================================================================*/     

void fill_w_tabulated(CPUComplex<CMPTYPE>* w_tabulated){
  int i,j;
  double x,y;
  ccomplex z;
  for(i=0;i<LENGTH;i++){
    y = WIDTH*(i-1);
    for(j=0;j<LENGTH;j++){
      x = WIDTH*(j-1);
      z = ccomplex(x,y);
      w_tabulated[i*LENGTH+j] = Faddeeva::w(z);
    }
  }
  return;
}

/*===============================================================================                   
 W_FUNCTION calculates the Faddeeva function, also known as the complex                
 probability integral, for complex arguments. For |z| < 6, it uses a six-point 
 interpolation scheme based on pre-tabulated data that is accurate to          
 O(10^-3). For |z| > 6, it uses a three-term asymptotic approximation that is                 
 accurate to O(10^-6).                           
===============================================================================*/ 
CPUComplex<CMPTYPE> w_function(CPUComplex<CMPTYPE> z, CPUComplex<CMPTYPE>* w_tabulated){
  CPUComplex<CMPTYPE> w;
  bool lower = (imag(z)<0);
  if(lower) z = Conjugate(z);
  //TODO: if(!lower){}else{...} should be faster but isn't
  if(abs(Norm(z)) < RADIUS){
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
///////////////////////////////////
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
////////////////////////////////////////////
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
  }
  return w;
}

