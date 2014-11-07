#include <optix_world.h>
#include "helpers.h"
#include "commonStructs.h"
#include <curand_kernel.h>
#include <math.h>

using namespace optix;

rtDeclareVariable(float,         scene_epsilon, , );
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

rtCallableProgram(double, xs_eval, (double));

#if defined(__MANY__)
__device__ unsigned long long intpow(int x, int n){
  unsigned long long y=1;
  for(int j=0;j<n;j++)
    y*=x; 
  return y;
}
#endif

RT_PROGRAM void generate_ray()
{
  unsigned nid = input_id_buffer[launch_index]%launch_dim;
  if(output_live_buffer[nid]){
  float phi = input_dir_a_buffer[nid];
  float mu  = input_dir_p_buffer[nid]; 
  float3 ray_origin = make_float3(input_pos_x_buffer[nid],input_pos_y_buffer[nid],input_pos_z_buffer[nid]);
  float3 ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 

  double xs,E;
  E = 2;
  xs = xs_eval(E); 
  printf("xs_eval(%g)=%g\n", E, xs);
#if defined(__MANY__)
  float3 direction_z = make_float3(0.f,0.f,1.f);
  optix::Ray rayz = optix::make_Ray(ray_origin, direction_z, only_one_ray_type, scene_epsilon, RT_DEFAULT_MAX);
#endif
  PerRayData_radiance prd;
  optix::Ray ray = optix::make_Ray(ray_origin, ray_direction, only_one_ray_type, scene_epsilon, RT_DEFAULT_MAX);
  prd.closestID = 0;
#if defined(__MANY__)
  prd.key = 0;
  prd.out = 0;
  prd.hits = 0;
  rtTrace(top_object, ray, prd);
  prd.closestID = prd.hitID;
  prd.closest_t = prd.t_hit;

  prd.key = 0;
  prd.out = 0;
  prd.hits = 0;
  rayz.origin = ray.origin + scene_epsilon*0.3*ray.direction;
  rtTrace(top_object, rayz, prd);
  rayz.origin = rayz.origin + prd.t_hit*rayz.direction;
  while(1-prd.out){
    rtTrace(top_object, rayz, prd);
    rayz.origin = rayz.origin + prd.t_hit*rayz.direction;
  }
  prd.current = prd.key/intpow(var_num,prd.hits-1);
#else
  prd.closest_t = MAX_LENGTH;
  prd.closestp_t = MAX_LENGTH;
  prd.closestID = 0;
  prd.closestpID = 0; 
  rtTrace(top_object, ray, prd);
  prd.current = prd.closestpID;
#endif
//#if defined(__PRINTTRACK__)
//  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
//         launch_index,prd.current,ray.origin.x,ray.origin.y,ray.origin.z);
//  ray.origin = ray.origin+(prd.closest_t+scene_epsilon*0.5)*ray.direction;
//  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
//         launch_index,1111,ray.origin.x,ray.origin.y,ray.origin.z); 
//#endif
  output_closest_buffer[nid] = prd.closest_t+scene_epsilon*0.5;
  output_current_buffer[nid] = prd.imat*(1-(0==prd.current));
  // if prd.current is found to be 0, the neutron leaks
  output_live_buffer[nid] = !(0==prd.current);
  //TODO: not determined whether closestID is needed
}
}


RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
  output_closest_buffer[launch_index] = 0.f;
}
