#include <optix_world.h>
#include "commonStructs.h"

using namespace optix;

rtDeclareVariable(float4,  sphere, , );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float,         scene_epsilon, , );
//rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
//rtDeclareVariable(int, geometryInstanceID, ,);
#if !defined(__MANY__)
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
#endif


RT_PROGRAM void sphere_intersect(int)
{
  float3 center = make_float3(sphere);
  float3 O = ray.origin - center;
  float3 D = ray.direction;
  float radius = sphere.w;

  float b = dot(O, D);
  float c = dot(O, O)-radius*radius;
  float disc = b*b-c;
  if(disc > 0.0f){
    disc = sqrtf(disc);
    c      = -b-disc; //t1
    radius = -b+disc; //t2
#if !defined(__MANY__)
    prd.onlyonce = (c*radius<=scene_epsilon);
#endif
    if( rtPotentialIntersection(c     /*t1*/)) {
      rtReportIntersection(0);
    } 
    if( rtPotentialIntersection(radius/*t2*/)) {
      rtReportIntersection(0);
    }
  }
}


RT_PROGRAM void sphere_bounds (int, float result[6])
{
  const float3 cen = make_float3( sphere );
  const float3 rad = make_float3( sphere.w );

  optix::Aabb* aabb = (optix::Aabb*)result;
  
  if( rad.x > 0.0f  && !isinf(rad.x) ) {
    aabb->m_min = cen - rad;
    aabb->m_max = cen + rad;
  } else {
    aabb->invalidate();
  }
}

