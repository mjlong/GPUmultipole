#include <optix_world.h>
#include "helpers.h"
#include "commonStructs.h"
#include <curand_kernel.h>
#include <math.h>
#include "global.h"

#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double
#define CMPTYPE2 double2
#endif

using namespace optix;

rtDeclareVariable(float,         var_R1, , );
rtDeclareVariable(float,         var_Hh, , );
rtDeclareVariable(unsigned,      var_num, , );

rtBuffer<unsigned, 1>           input_id_buffer;

rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<float, 1>              input_pos_y_buffer;
rtBuffer<float, 1>              input_pos_z_buffer;

rtBuffer<float, 1>              input_dir_p_buffer;
rtBuffer<float, 1>              input_dir_a_buffer;

rtBuffer<float,    1> output_closest_buffer;
rtBuffer<unsigned, 1> output_current_buffer;
rtBuffer<unsigned, 1> output_live_buffer;

rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(unsigned int,  only_one_ray_type, , );

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

rtCallableProgram(void, xs_eval, (int, CMPTYPE, CMPTYPE, CMPTYPE*,CMPTYPE*,CMPTYPE* ));
rtCallableProgram(void, locate,  (float3, float3, float*, unsigned*, unsigned* ));

RT_PROGRAM void generate_ray()
{
  unsigned nid = input_id_buffer[launch_index]%launch_dim;
  if(output_live_buffer[nid]){
  float phi = input_dir_a_buffer[nid];
  float mu  = input_dir_p_buffer[nid]; 
  float3 ray_origin = make_float3(input_pos_x_buffer[nid],input_pos_y_buffer[nid],input_pos_z_buffer[nid]);
  float3 ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 
  float d;
  unsigned icell, imat;

  double E=2;
  double sigT,sigA,sigF;
  xs_eval(0,E,sqrt(300.*KB),&sigT,&sigA,&sigF); 
  printf("xs_eval(%g)=%g,%g,%g\n",E,sigT,sigA,sigF);

  locate(ray_origin,ray_direction, &d, &imat, &icell);
#if defined(__PRINTTRACK__)
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,icell,ray_origin.x,ray_origin.y,ray_origin.z);
  ray_origin = ray_origin+d*ray_direction;
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,1111,ray_origin.x,ray_origin.y,ray_origin.z); 
#endif


  output_closest_buffer[nid] = d;
  output_current_buffer[nid] = imat*(1-(0==icell));
  output_live_buffer[nid] = !(0==icell);
}
}

RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
  output_closest_buffer[launch_index] = 0.f;
}
