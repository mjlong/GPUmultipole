#include <optix_world.h>
#include "commonStructs.h"
#include "helpers.h"

struct PerRayData_radiance
{
  float3 result;
  float importance;
  int depth;
};

struct PerRayData_shadow
{
  float3 attenuation;
};


rtDeclareVariable(unsigned int,      radiance_ray_type, , );
rtDeclareVariable(unsigned int,      shadow_ray_type, , );
rtDeclareVariable(rtObject,          top_object, , );
rtDeclareVariable(rtObject,          top_shadower, , );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
rtDeclareVariable(PerRayData_shadow,   prd_shadow, rtPayload, );

/*
static __device__ void phongShadowed()
{
  // this material is opaque, so it fully attenuates all shadow rays
  prd_shadow.attenuation = optix::make_float3(0);
  rtTerminateRay();
}
*/

static
__device__ void phongShade( float3 p_normal)
{
  float3 hit_point = ray.origin + t_hit * ray.direction;
  printf("ray from (%g,%g,%g) travels along (%g,%g,%g) for %g to (%g,%g,%g)\n", 
         ray.origin.x, ray.origin.y, ray.origin.z, ray.direction.x, ray.direction.y, ray.direction.z,t_hit, 
         hit_point.x, hit_point.y, hit_point.z);  

  float3 result = make_float3(0.f,0.f,0.f);
  prd.result = result;
}

