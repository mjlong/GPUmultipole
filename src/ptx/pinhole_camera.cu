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

rtBuffer<CMPTYPE, 1>            input_energy_buffer;
rtBuffer<float, 1>              input_pos_x_buffer;
rtBuffer<float, 1>              input_pos_y_buffer;
rtBuffer<float, 1>              input_pos_z_buffer;

rtBuffer<float, 1>              input_dir_p_buffer;
rtBuffer<float, 1>              input_dir_a_buffer;

rtBuffer<unsigned, 1>           output_live_buffer;
rtBuffer<unsigned, 1>           output_terminated_buffer;
rtBuffer<unsigned, 1>           output_spectrum_buffer;


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


__device__ void neutron_sample(unsigned* live, CMPTYPE* energy, float3* origin, float* mu, float* phi, curandState* localstate){
  *live = 1u;
  *energy = STARTENE;
  *phi =   2*PI*curand_uniform(localstate);
  *mu  = -1.f+2*curand_uniform(localstate); 
  *origin =  make_float3(0.5f+0.00*curand_uniform(localstate),
                                     0.5f+0.00*curand_uniform(localstate),
                                     0.5f+0.00*curand_uniform(localstate));
  output_terminated_buffer[launch_index]+=1;
}

RT_PROGRAM void generate_ray()
{
  unsigned icell, imat,isotope,live,stay;
#if defined(__PRINTTRACK__)
  int nid = output_live_buffer[launch_index]/2;
  live = output_live_buffer[launch_index]%2;
#else
  live = output_live_buffer[launch_index];
#endif
  float d;
  float mu,phi;
  CMPTYPE localenergy;
  double sigT,sigA,sigF,
         sigTsum,sigAsum,sigFsum;
  float3 ray_origin, ray_direction;
  curandState localstate = input_random[launch_index];
  localenergy = input_energy_buffer[launch_index];
  ray_origin.x = input_pos_x_buffer[launch_index];
  ray_origin.y = input_pos_y_buffer[launch_index];
  ray_origin.z = input_pos_z_buffer[launch_index];
  mu  = input_dir_p_buffer[launch_index]; 
  phi = input_dir_a_buffer[launch_index];
  ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 
  //stay = 0;
//loop over GPU steps
  for(unsigned istep=0; istep<devstep; istep++){
  //if(!stay){
    locate(ray_origin,ray_direction, &d, &imat, &icell);
  //  printf("rt invoked\n");
  //}
  imat = imat*(1-(0==icell));
  live = !(0==icell);

  if(!live){
#if defined(__PRINTTRACK__)
    if(__PRINTTRACK__){
    printf("%7d,%3d,%+.7e, %+.7e, %+.7e, %.14e leaked\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy); 
    }
    nid += launch_dim;
#endif
    neutron_sample(&live, &localenergy, &ray_origin, &mu, &phi, &localstate);
    ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 
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
  float s = -log(curand_uniform(&localstate))/1E8;//1.4;
  //invoking xs_eval not only messes nid but makes phi unchanged
  //printf("%7d,s=%g,mu=%g, phi=%g\n",nid,s,mu,phi);
  //stay = (d>s);
  s = (d<s)*d + (d>=s)*s;
    
//update tally
  output_spectrum_buffer[search_bin(localenergy)*launch_dim+launch_index]+=1; 
  live = (localenergy > ENDENERG);

  //localenergy = localenergy*live + STARTENE*(1u-live);
  if(live){
    ray_origin = ray_origin + s*ray_direction;
  }
  else{
    //stay = 0;
#if defined(__PRINTTRACK__)
    if(__PRINTTRACK__){
    printf("%7d,%3d,%+.7e, %+.7e, %+.7e, %.14e stopped\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy); 
    }
    nid += launch_dim;
#endif
    neutron_sample(&live, &localenergy, &ray_origin, &mu, &phi, &localstate); 
    ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 
  }
  }//end for istep
  input_random[launch_index]  = localstate;
#if defined(__PRINTTRACK__)
  output_live_buffer[launch_index] = 2*nid+live;
#else
  output_live_buffer[launch_index] = live;
#endif
  input_pos_x_buffer[launch_index] = ray_origin.x;
  input_pos_y_buffer[launch_index] = ray_origin.y;
  input_pos_z_buffer[launch_index] = ray_origin.z;
  input_dir_p_buffer[launch_index] = mu;
  input_dir_a_buffer[launch_index] = phi;
  input_energy_buffer[launch_index] = localenergy;
}

RT_PROGRAM void remaining_ray()
{
  unsigned icell, imat,isotope,live,stay;
#if defined(__PRINTTRACK__)
  int nid = output_live_buffer[launch_index]/2;
  live = output_live_buffer[launch_index]%2;
// terminated neutrons in remaining batch is not counted by 
// DeviceMem.grid_terminated_neutrons[launch_dim] but 
// count how many are live before remaining() batch
  output_live_buffer[launch_index] = live; 
#else
  live = output_live_buffer[launch_index];
#endif
  float d;
  float mu,phi;
  CMPTYPE localenergy;
  double sigT,sigA,sigF,
         sigTsum,sigAsum,sigFsum;
  float3 ray_origin, ray_direction;
  curandState localstate = input_random[launch_index];
  localenergy = input_energy_buffer[launch_index];
  ray_origin.x = input_pos_x_buffer[launch_index];
  ray_origin.y = input_pos_y_buffer[launch_index];
  ray_origin.z = input_pos_z_buffer[launch_index];
  mu  = input_dir_p_buffer[launch_index]; 
  phi = input_dir_a_buffer[launch_index];
  ray_direction = make_float3(sqrt(1.f-mu*mu)*cos(phi),sqrt(1.f-mu*mu)*sin(phi),mu); 
  //stay = 0;
//loop over remaining neutron life 
  while(live){
  //if(!stay){
    locate(ray_origin,ray_direction, &d, &imat, &icell);
  //  printf("rt invoked\n");
  //}
  imat = imat*(1-(0==icell));
  live = !(0==icell);

  if(!live){
#if defined(__PRINTTRACK__)
    if(__PRINTTRACK__){
    printf("[r]%4d,%3d,%+.7e, %+.7e, %+.7e, %.14e leaked\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy); 
    }
    nid += launch_dim;
#endif
    break;
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
    printf("[r]%4d,%3d,%+.7e, %+.7e, %+.7e, %.14e %.14e %.14e %.14e\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy, sigTsum,sigAsum,sigFsum); 
  }
#endif
  localenergy = localenergy*curand_uniform(&localstate);
  float s = -log(curand_uniform(&localstate))/1E8;//1.4;
  //stay = (d>s);
  s = (d<s)*d + (d>=s)*s;
  //update tally
  output_spectrum_buffer[search_bin(localenergy)*launch_dim+launch_index]+=1; 
  live = (localenergy > ENDENERG);
  //localenergy = localenergy*live + STARTENE*(1u-live);
  if(live){
    ray_origin = ray_origin + s*ray_direction;
  }
  else{
#if defined(__PRINTTRACK__)
    if(__PRINTTRACK__){
    printf("[r]%4d,%3d,%+.7e, %+.7e, %+.7e, %.14e stopped\n",
            nid, imat,
            ray_origin.x, ray_origin.y, ray_origin.z,
            localenergy); 
    }
    nid += launch_dim;
#endif
    break;
  }
  }//end while(live) 

}

RT_PROGRAM void exception()
{
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception 0x%X at launch index (%d)\n", code, launch_index );
}
