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

rtDeclareVariable(unsigned,      devstep, , );

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

__device__ __constant__ CMPTYPE spectrumbins[]={
  20000000.0,
  1000000.0,
  100000.0,
  10000.0,
  1000.0,
  100.0,
  50.0,
  25.0,
  10.0,
  6.0,
  4.0,
  1.0,
  0.625,
  0.3,
  0.1,
  0.03,
  0.00001
};

__device__ unsigned search_bin(CMPTYPE energy){
  for(int i=0;i<NUM_BINS;i++){
    if( (spectrumbins[i]>=energy)&&(spectrumbins[i+1]<energy) ) 
      return i;
  }
  return 0;
}


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
  unsigned icell, imat,isotope,live;
#if defined(__PRINTTRACK__)
  int nid = launch_index;
#endif
  float d;

  CMPTYPE localenergy;
  double sigT,sigA,sigF,
         sigTsum,sigAsum,sigFsum;
  float3 ray_origin, ray_direction;
  curandState localstate = input_random[launch_index];

  neutron_sample(&live, &localenergy, &ray_origin, &ray_direction, &localstate);

//loop over GPU steps
  for(unsigned istep=0; istep<devstep; istep++){
  locate(ray_origin,ray_direction, &d, &imat, &icell);
  imat = imat*(1-(0==icell));
  live = !(0==icell);

  if(!live){
    neutron_sample(&live, &localenergy, &ray_origin, &ray_direction, &localstate);
#if defined(__PRINTTRACK__)
    printf("leaked\n");
    nid += launch_dim;
#endif
  }
//
//Evaluate cross section and print (id,imat,position,E,sigT,sigA,sigF
//
  sigTsum = 0;
  sigAsum = 0;
  sigFsum = 0;
  for(isotope=mat_offsets[imat];isotope<mat_offsets[imat+1];isotope++ ){
    //xs_eval(mat_isotopes[isotope],localenergy,sqrt(300.*KB),&sigT,&sigA,&sigF); 
    sigTsum += sigT*mat_densities[isotope];
    sigAsum += sigA*mat_densities[isotope];
    sigFsum += sigF*mat_densities[isotope];
  }
#if defined(__PRINTTRACK__)
  if(__PRINTTRACK__){
    printf("%7d,%3d,%+.7e, %+.7e, %+.7e, %.14e %.14e %.14e %.14e\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy, sigTsum,sigAsum,sigFsum); 
  }
#endif
  localenergy = localenergy*curand_uniform(&localstate);
  float s = -log(curand_uniform(&localstate))/1.4;
  s = (d<s)*d + (d>=s)*s;
//update tally
  int iE;
  iE = search_bin(localenergy); 
  live = (localenergy > ENDENERG);
  //localenergy = localenergy*live + STARTENE*(1u-live);
  if(live){
    ray_origin = ray_origin + s*ray_direction;
  }
  else{
    printf("stopped\n");
    neutron_sample(&live, &localenergy, &ray_origin, &ray_direction, &localstate); 
#if defined(__PRINTTRACK__)
    nid += launch_dim;
#endif
  }

  }//end for istep
}

RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
}
