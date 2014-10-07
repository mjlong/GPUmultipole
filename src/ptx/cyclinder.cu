#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>
#include "commonStructs.h"

using namespace optix;

rtDeclareVariable(float3, p1, , );
rtDeclareVariable(float3, p2, , );
rtDeclareVariable(float , R , , );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float,         scene_epsilon, , );
rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(int, geometryInstanceID, ,);
#if !defined(__MANY__)
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
#endif

RT_PROGRAM void cylinder_intersect(int)
{
//  if(3==launch_index && 23==geometryInstanceID)
//    printf("x=%g,y=%g\n",p1.x,p1.y);
  float t1=0.f;
  float t2=0.f;
  float t3=0.f;
  float t4=0.f;
  float tmin=0.f;
  float tmax=0.f;
  float3 m = normalize(p2-p1);
  float md = dot(m,ray.direction);
  float mrp= dot(m,ray.origin-p1);
  float A = 1.f-md*md;
  float B = dot(ray.origin-p1,ray.direction) -md*mrp;
  float C = dot(ray.origin-p1,ray.origin-p1) - mrp*mrp - R*R;
  tmin = B*B-A*C;
//  if(3==launch_index && 23==geometryInstanceID)
//    printf("A=%g,B=%g,C=%g,delta=%g\n",A,B,C,tmin);

  if(0==A){
    //t3 = length(mrp*m+p1-r0)-R;
    //if(t3>0){
    //  t1 = 1.f/0.f;
    //  t2 = 1.f/0.f;
    //} else{
    //  t1 = -1.f/0.f;
    //  t2 =  1.f/0.f;
    //}
    t1 = (((length(mrp*m+p1-ray.origin)>R)<<1) -1)*1.f/A;
    t2 = 1.f/A;
  } else if(tmin>0){
    t1 = (-B-sqrtf(B*B-A*C))/A;
    t2 = (-B+sqrtf(B*B-A*C))/A;
  } else{
    t1 = 1.f/0.f;
    t2 = 1.f/0.f;
  }
  A = -mrp/md;
  B = dot(p2-ray.origin,m)/md;
  t3 = min(A,B); 
  t4 = max(A,B);
  tmin = max(t1,t3);
  tmax = min(t2,t4);
  //printf("i'm thread %d in intersection with %d,t1=%g,t2=%g,t3=%g,t4=%g, tmin=%g, tmax=%g\n", launch_index,geometryInstanceID,t1,t2,t3,t4,tmin,tmax);
  if(tmin <= tmax) {
//  if(3==launch_index && 23==geometryInstanceID)
//    printf("t1=%g,t2=%g,t3=%g,t4=%g,tmin=%g,tmax=%g\n",t1,t2,t3,t4,tmin,tmax);
#if !defined(__MANY__)
    prd.onlyonce = (1==((tmax<scene_epsilon) + (tmin<scene_epsilon))); 
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

RT_PROGRAM void cylinder_bounds (int, float result[6])
{
  optix::Aabb* aabb = (optix::Aabb*)result;
  float3 m = normalize(p2-p1);
  float3 parallel = normalize(make_float3(1.f,1.f,0.f-(m.x+m.y)/m.z));
  aabb->set(p1-R*parallel, p2+R*parallel);
}
