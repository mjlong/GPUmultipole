#include "multipole.h"
#include "CComplex.h"
#include <optix.h>

#if defined (__QUICKWC) 
extern __constant__ CMPTYPE2 constwtable[LENGTH*LENGTH];
#endif

#if defined (__MITW)
#include "Faddeeva.h"
#define w_function Faddeeva::w
#endif
rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<unsigned int, 1>       dev_integers;
rtBuffer<CComplex<double>, 1>   mpdata;

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

RT_CALLABLE_PROGRAM double xs_eval(double E){
  CComplex<double> z=mpdata[launch_index];
  printf("z=%g i%+g\n",real(z),imag(z));
  //CComplex<double>z1 = CComplex<double>(constwtable[0].x,constwtable[0].y);
  CComplex<double> z1 = w_function(z,0.00001);
  printf("z1=%g i%+g\n",real(z1),imag(z1));

  return 1/E;
}

