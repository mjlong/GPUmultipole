//#include "multipole.h"
#include "CComplex.h"
#include <optix.h>

#if defined (__QUICKWC) 
extern __constant__ CMPTYPE2 constwtable[];
#endif

#if defined (__MITW)
#include "Faddeeva.h"
#define w_function Faddeeva::w
#endif
rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<unsigned int, 1>       dev_integers;
rtBuffer<CComplex<double>, 1>   mpdata;
rtBuffer<CComplex<double>, 1>   wtable_buffer;

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

RT_CALLABLE_PROGRAM double xs_eval(double E){
  CComplex<double> z=mpdata[launch_index];
  printf("z=%g%+gi\n",real(z),imag(z));
#if defined(__QUICKWC)
  CComplex<double>z1 = CComplex<double>(constwtable[34].x,constwtable[34].y);
  CComplex<double>z2 = wtable_buffer[34];
  printf("constant in cu,       z1=%g%+gi\n",real(z1),imag(z1));
  printf("from global,          z2=%g%+gi\n",real(z2),imag(z2));
#endif
#if defined(__MITW)
  CComplex<double> zzzz = CComplex<double>(2.2,1.3);
  printf("zzzz=%g%+gi\n",real(zzzz),imag(zzzz));

  //printf("fakew(%g)=%g\n",3.2, fakew(3.2));
  //double xx = fakew(2.2);
  CComplex<double> z1 = w_function(z,0.00001);
  //printf("F(%g%+gi)=%g%+gi\n",real(z),imag(z),real(z1),imag(z1));
#endif
  return 1/E;
}

