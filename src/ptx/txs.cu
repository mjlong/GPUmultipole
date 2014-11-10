#include "multipole.h"
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
  printf("z=%g i%+g\n",real(z),imag(z));
#if defined(__QUICKWC)
  CComplex<double>z1 = CComplex<double>(constwtable[34].x,constwtable[34].y);
  CComplex<double>z2 = wtable_buffer[34];
#endif
#if defined(__MITW)
  CComplex<double> z1 = w_function(z,0.00001);
#endif
  printf("z1=%g i%+g\n",real(z1),imag(z1));
  printf("z2=%g i%+g\n",real(z2),imag(z2));

  return 1/E;
}

