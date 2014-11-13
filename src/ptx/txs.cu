//#include "multipole.h"
#include "CComplex.h"
#include <optix.h>

#if defined(__MITW)
#include "Faddeeva.h"
#define w_function Faddeeva::w
#endif

#if defined(__QUICKW)
#include "QuickW.h"
rtBuffer<CComplex<double>, 1>   wtable_buffer;
#endif

#if defined(__FOURIERW)
#include "fourierw.h"
#endif


rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<unsigned int, 1>       dev_integers;
rtBuffer<CComplex<double>, 1>   mpdata;

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

RT_CALLABLE_PROGRAM double xs_eval(double E){
#if defined(__MITW)
  
#endif

  return 1/E;
}

