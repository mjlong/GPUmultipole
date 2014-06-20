#include "QuickW.h"

double b = 0.275255128608410950901357962647054304017026259671664935783653;
double d = 2.724744871391589049098642037352945695982973740328335064216346;
double a = 0.512424224754768462984202823134979415014943561548661637413182;
double c = 0.051765358792987823963876628425793170829107067780337219430904;

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

__device__ void fill_w_tabulated(CComplex* w_tabulated, unsigned id){
  double x,y;
  CComplex z;
  y = WIDTH*(id/LENGTH-1);
  x = WIDTH*(id%LENGTH-1);
  z = CComplex(x,y);
  w_tabulated[id] = Faddeeva::w(z);
  return;
}


/*===============================================================================                   
 W_FUNCTION calculates the Faddeeva function, also known as the complex                
 probability integral, for complex arguments. For |z| < 6, it uses a six-point 
 interpolation scheme based on pre-tabulated data that is accurate to          
 O(10^-3). For |z| > 6, it uses a three-term asymptotic approximation that is                 
 accurate to O(10^-6).                           
===============================================================================*/ 
__device__ CComplex w_function(CComplex z, CComplex* w_tabulated){
  double  p;           // interpolation factor on real axis                                   
  double  q;           // interpolation factor on imaginary axis                                  
  double  pp, qq, pq;  // products of p and q                                         
  double  a_l;         // coefficient for left point                                   
  double  a_c;         // coefficient for center point                                         
  double  a_b;         // coefficient for bottom point    
  double  a_r;         // coefficient for right point                                             
  double  a_t;         // coefficient for top point  

  int l;               //interpolation index for real axis
  int m;               //interpolation index for imaginary axis
  
  CComplex w;
  
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

    //Coefficients for interpolation
    a_b = 0.5*(qq - q);         //bottom
    a_l = 0.5*(pp - p);         //left
    a_c = 1.0 + pq - pp - qq;   //center
    a_r = 0.5*(pp + p) - pq;    //right
    a_t = 0.5*(qq + q) - pq;    //top
  
    // Use six-point interpolation to calculate real and imaginary parts
    l++;
    m++;
    w =  
      a_b*w_tabulated[(m-1)*LENGTH+l] + 
      a_l*w_tabulated[m*LENGTH + l-1] +
      a_c*w_tabulated[m*LENGTH + l  ] +
      a_r*w_tabulated[m*LENGTH + l+1] +
      a_t*w_tabulated[(m+1)*LENGTH+l] +
      pq *w_tabulated[(m+1)*LENGTH+l+1];
    if(real(z)<0) 
      w = Conjugate(w);
  }
  else
    w = ONEI * z * (a/(z*z - b) + c/(z*z - d));

  return w;
}
