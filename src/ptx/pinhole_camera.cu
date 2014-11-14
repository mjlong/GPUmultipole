#include <optix_world.h>
#include "helpers.h"
#include "commonStructs.h"
#include <curand_kernel.h>
#include <math.h>

#include "simulation.h"
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

rtBuffer<unsigned, 1>           mat_offsets;
rtBuffer<unsigned, 1>           mat_isotopes;
rtBuffer<float, 1>              mat_densities;

rtBuffer<curandState, 1>           input_random;

rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(unsigned int,  only_one_ray_type, , );

rtDeclareVariable(unsigned int, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, launch_dim,   rtLaunchDim, );

rtCallableProgram(void, xs_eval, (int, CMPTYPE, CMPTYPE, CMPTYPE*,CMPTYPE*,CMPTYPE* ));
rtCallableProgram(void, locate,  (float3, float3, float*, unsigned*, unsigned* ));

__device__ void neutron_sample(unsigned* live, CMPTYPE* energy, float3* origin, float3* direction, curandState* localstate){
  *live = 1u;
  *energy = STARTENE;
  float phi =   2*PI*curand_uniform(localstate);
  float mu  = -1.f+2*curand_uniform(localstate); 
  *origin =  make_float3(0.5f+0.00*curand_uniform(localstate),
                                     0.5f+0.00*curand_uniform(localstate),
                                     0.5f+0.00*curand_uniform(localstate));
  *direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu);
}


RT_PROGRAM void generate_ray()
{
  curandState localstate = input_random[launch_index];  
  float3 ray_origin, ray_direction;
  CMPTYPE localenergy;
  float d;
  unsigned icell, imat,isotope,live;

  neutron_sample(&live, &localenergy, &ray_origin, &ray_direction, &localstate);

  locate(ray_origin,ray_direction, &d, &imat, &icell);
  imat = imat*(1-(0==icell));
  live = !(0==icell);

  if(!live)
    neutron_sample(&live, &localenergy, &ray_origin, &ray_direction, &localstate);

  double sigT,sigA,sigF,
         sigTsum,sigAsum,sigFsum;

  sigTsum = 0;
  sigAsum = 0;
  sigFsum = 0;
  for(isotope=mat_offsets[imat];isotope<mat_offsets[imat+1];isotope++ ){
    xs_eval(mat_isotopes[isotope],localenergy,sqrt(300.*KB),&sigT,&sigA,&sigF); 
    sigTsum += sigT*mat_densities[isotope];
    sigAsum += sigA*mat_densities[isotope];
    sigFsum += sigF*mat_densities[isotope];
  }

  printf("xs_eval(%g)=%g,%g,%g\n",localenergy,sigTsum,sigAsum,sigFsum);

#if defined(__PRINTTRACK__)
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,icell,ray_origin.x,ray_origin.y,ray_origin.z);
  ray_origin = ray_origin+d*ray_direction;
  printf("%3d, %3d, %+18.12e,%+18.12e,%+18.12e\n",
         launch_index,1111,ray_origin.x,ray_origin.y,ray_origin.z); 
#endif


}

RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
}
