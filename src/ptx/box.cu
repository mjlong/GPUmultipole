#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>
#include "commonStructs.h"

using namespace optix;

rtDeclareVariable(float3, boxmin, , );
rtDeclareVariable(float3, boxmax, , );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float,         scene_epsilon, , );
//rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
//rtDeclareVariable(int, geometryInstanceID, ,);
#if !defined(__MANY__)
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
#endif


RT_PROGRAM void box_intersect(int)
{
  float3 t0 = (boxmin - ray.origin)/ray.direction;
  float3 t1 = (boxmax - ray.origin)/ray.direction;
  float3 near = fminf(t0, t1);
  float3 far = fmaxf(t0, t1);
  float tmin = fmaxf( near );
  float tmax = fminf( far );
  //
  //printf("i'm thread %d in intersection, tmin=%g, tmax=%g\n", launch_index,tmin,tmax);
  if(tmin <= tmax) {
#if !defined(__MANY__)
    prd.onlyonce = (tmax*tmin<=scene_epsilon); 
#endif
    if( rtPotentialIntersection( tmax ) ) {
      //printf("tmax passed. i'm thread (%d,%d), ray from (%g,%g,%g) travels %g along (%g,%g,%g)\n", 
      //       threadIdx.x+blockDim.x*blockIdx.x, launch_index, ray.origin.x, ray.origin.y, ray.origin.z, 
      //       tmax*length(ray.direction), ray.direction.x, ray.direction.y, ray.direction.z);
      rtReportIntersection(0);
    }
    if( rtPotentialIntersection( tmin ) ) {
      //printf("tmin passed. i'm thread (%d,%d), ray from (%g,%g,%g) travels %g along (%g,%g,%g)\n", 
      //        threadIdx.x+blockDim.x*blockIdx.x, launch_index, ray.origin.x, ray.origin.y, ray.origin.z, 
      //        tmin*length(ray.direction), ray.direction.x, ray.direction.y, ray.direction.z);
      rtReportIntersection(0);
      //if(rtReportIntersection(0))
      //  check_second = true;
   } 
 }
}

RT_PROGRAM void box_bounds (int, float result[6])
{
  optix::Aabb* aabb = (optix::Aabb*)result;
  aabb->set(boxmin, boxmax);
}
