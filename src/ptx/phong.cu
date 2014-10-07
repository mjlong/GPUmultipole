#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include "commonStructs.h"

using namespace optix;

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(unsigned,var_num, ,);

//rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(int, geometryInstanceID, , );
rtDeclareVariable(PerRayData_radiance, prd,rtPayload , );

#if defined(__MANY__)
RT_PROGRAM void closest_hit_radiance()
{
  prd.hitID = geometryInstanceID;
  prd.t_hit  = t_hit;
//  if(2==launch_index&&1.f==ray.direction.z)
//    printf("%d hitting %d, t=%g\n", launch_index, geometryInstanceID, t_hit);
  if(prd.key%var_num==geometryInstanceID){
    prd.key = (prd.key-geometryInstanceID)/var_num;
    prd.hits--;
  }
  else{
    prd.key = prd.key*var_num+geometryInstanceID; 
    prd.hits++;
  }
}
#else
RT_PROGRAM void any_hit_shadow()
{
  prd.current = geometryInstanceID; 
  if(prd.closest_t> t_hit){
    prd.closest_t = t_hit;
    prd.closestID = geometryInstanceID;
  }
  if(prd.onlyonce && (prd.closestp_t > t_hit)){
    prd.closestp_t = t_hit;
    prd.closestpID = geometryInstanceID;
  }
  rtIgnoreIntersection();
}
#endif
