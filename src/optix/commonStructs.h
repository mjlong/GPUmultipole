#pragma once
#include <optixu/optixu_vector_types.h>
#define MAX_LENGTH 9999.9
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#if defined(__BVH__)
#define BUILDER "Trbvh"
#define TRAVERSER "Bvh"
#else
#define BUILDER "NoAccel"
#define TRAVERSER "NoAccel"
#endif

struct PerRayData_radiance
{
  unsigned closestID;
  float    closest_t;
  unsigned current;
  unsigned imat;
#if defined(__MANY__)
  unsigned long long key;
  unsigned out;
  unsigned hits;
  unsigned hitID;
  float    t_hit;
#else
  unsigned onlyonce; 
  unsigned closestpID;
  float    closestp_t;
#endif
};
