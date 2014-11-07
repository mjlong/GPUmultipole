#include <optix.h>
rtBuffer<float, 1>              input_pos_x_buffer;
rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

using namespace optix;

RT_CALLABLE_PROGRAM double xs_eval(double E){
  input_pos_x_buffer[3-launch_index]=0.3;
  printf("x[%d]=%g\n",launch_index, input_pos_x_buffer[launch_index]);
  return 1/E;
}

