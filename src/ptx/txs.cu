#include "multipole.h"
#include <optix.h>
rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<unsigned int, 1>       dev_integers;

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

RT_CALLABLE_PROGRAM double xs_eval(double E){
  return 1/E;
}

