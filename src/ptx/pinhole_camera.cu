#include <optix_world.h>
#include "helpers.h"
#include "commonStructs.h"
#include <curand_kernel.h>
#include <math.h>

using namespace optix;
#define NUM_COL 64 

rtDeclareVariable(float,         scene_epsilon, , );
rtDeclareVariable(float,         var_R1, , );
rtDeclareVariable(float,         var_Hh, , );
rtDeclareVariable(unsigned,      var_num, , );

rtBuffer<float3, 1>              input_pos_buffer;
rtBuffer<float3, 1>              input_dir_buffer;
rtBuffer<float, 1>              input_test_buffer;

rtBuffer<float, 1>    output_closest_buffer;
rtBuffer<unsigned, 1> output_current_buffer;
rtBuffer<unsigned, 1> output_next_buffer;
rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(unsigned int,  only_one_ray_type, , );

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

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
  curandState localstate;
  curand_init(launch_index,0,0,&localstate);
  printf("test =%f\n",input_test_buffer[launch_index]);
  int iCol = 0;
  float lambda = 0.f;
#if defined(__HEXPRISM__)
  float phi = 2*PI*curand_uniform(&localstate);
  float mu  = -1.f+2*curand_uniform(&localstate);
  //float3 ray_origin = make_float3(0.48f,0.f,0.f);//input_pos_buffer[launch_index];
  float3 ray_origin = make_float3(var_R1*cos(phi), var_R1*sin(phi), -var_Hh+2*var_Hh*curand_uniform(&localstate));
  phi = 2*PI*curand_uniform(&localstate);
  //float3 ray_direction = input_dir_buffer[launch_index];
  float3 ray_direction = normalize(make_float3(sqrt(1.f-mu)*cos(phi),sqrt(1.f-mu)*sin(phi),mu)); 
#else
  float phi = 2*PI*curand_uniform(&localstate);
  float mu  = -1.f+2*curand_uniform(&localstate);
  //float3 ray_origin = make_float3(0.48f,0.f,0.f);//input_pos_buffer[launch_index];
  float3 ray_origin = make_float3(var_R1*(0.5*cos(phi)+sqrt(0.5)/0.9), var_R1*(0.5*sin(phi)+sqrt(0.5)/0.9), -var_Hh+2*var_Hh*curand_uniform(&localstate));
  phi = 2*PI*curand_uniform(&localstate);
  //float3 ray_direction = input_dir_buffer[launch_index];
  float3 ray_direction = normalize(make_float3(sqrt(1.f-mu)*cos(phi),sqrt(1.f-mu)*sin(phi),mu)); 
#endif

#if defined(__MANY__)
  float3 direction_z = make_float3(0.f,0.f,1.f);
  optix::Ray rayz = optix::make_Ray(ray_origin, direction_z, only_one_ray_type, scene_epsilon, RT_DEFAULT_MAX);
#endif
  PerRayData_radiance prd;
  optix::Ray ray = optix::make_Ray(ray_origin, ray_direction, only_one_ray_type, scene_epsilon, RT_DEFAULT_MAX);
//start loop over iCollision
  //for(iCol=0;iCol<NUM_COL;iCol++){
  prd.closestID = 0;
  while((prd.closestID!=var_num-1 || lambda<prd.closest_t) && iCol<NUM_COL){
  lambda = max(0.48,-logf(curand_uniform(&localstate))/3.45*10.0);
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
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,/*iCol*/prd.current,ray.origin.x,ray.origin.y,ray.origin.z);
  //printf("%3d, %dth collision, history %llu, hits %3d, closest hit %3d,%+10.4e; in %3d\n",launch_index,iCol, prd.key,prd.hits,prd.closestID,prd.closest_t, prd.current);
  lambda = prd.closest_t;//lambda > prd.closest_t ? prd.closest_t : lambda;
  ray.origin = ray.origin + (lambda+scene_epsilon*0.5)*ray.direction;
#else

  prd.closest_t = MAX_LENGTH;
  prd.closestp_t = MAX_LENGTH;
  prd.closestID = 0;
  prd.closestpID = 0; 
  rtTrace(top_object, ray, prd);
  prd.current = prd.closestpID;
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,/*iCol*/prd.current,ray.origin.x,ray.origin.y,ray.origin.z);
  //printf("%3d, %dth collision, closest %3d,%+10.4e; closest' %3d,%+10.4e ,in %3d\n",
  //       launch_index,iCol, prd.closestID,prd.closest_t, prd.closestpID,prd.closestp_t, prd.current);
  lambda = prd.closest_t;//lambda > prd.closest_t ? prd.closest_t : lambda;
  ray.origin = ray.origin+(lambda+scene_epsilon*0.5)*ray.direction;
#endif
//  if(lambda < prd.closest_t){
//    //printf("%3d, %dth collision, remains in the cell.\n",launch_index, iCol);
//    phi = 2*PI*curand_uniform(&localstate);
//    mu  = -1.f+2*curand_uniform(&localstate);
//    ray.direction = normalize(make_float3(sqrt(1.f-mu)*cos(phi),sqrt(1.f-mu)*sin(phi),mu)); 
//  } else {
//    //printf("%3d, %dth collision, leaves current cell.\n",launch_index,iCol);
//  } 
  //printf("%3d, %dth collision, closestID = %3d,var_num-1=%3d\n",
  //       launch_index,iCol, prd.closestID,var_num-1);
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,iCol,ray.origin.x,ray.origin.y,ray.origin.z);
  iCol++; 
  }
//////end loop over iCol
  output_closest_buffer[launch_index] = prd.closest_t;
  output_current_buffer[launch_index] = prd.current;
}


RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
  output_closest_buffer[launch_index] = 0.f;
}
