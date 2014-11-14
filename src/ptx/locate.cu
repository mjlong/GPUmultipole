#include <optix_world.h>
#include "helpers.h"
#include "commonStructs.h"
#include <curand_kernel.h>
#include <math.h>

using namespace optix;

rtDeclareVariable(float,         scene_epsilon, , );

rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(unsigned int,  only_one_ray_type, , );


#if defined(__MANY__)
__device__ unsigned long long intpow(int x, int n){
  unsigned long long y=1;
  for(int j=0;j<n;j++)
    y*=x; 
  return y;
}
#endif

RT_CALLABLE_PROGRAM void locate(float3 ray_origin,float3 ray_direction, 
                             float* distance, unsigned* imat, unsigned* icell)
{
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
  *distance = prd.closest_t+scene_epsilon*0.5;
  *imat = prd.imat;
  *icell = prd.current;
}


