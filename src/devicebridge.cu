#include "CPUComplex.h"
#include "simulation.h"
#include "manmemory.h"

#include "devicebridge.h"
/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/
void print_results(unsigned num_src, unsigned num_bin,  MemStruct HostMem, float timems){
  
/*print collision cnt and time*/
  unsigned sum=0;
  for(int j=0;j<num_bin;j++){ 
    sum+=HostMem.spectrum[j];
    printf("%4d \n",HostMem.spectrum[j]);
  }
  printf("%u\n",HostMem.num_terminated_neutrons);
  printf("time elapsed:%g mus\n", timems*1000/sum);
  
  FILE *fp=NULL;
  fp = fopen("timelog","a+");
  fprintf(fp,"%-.6f,%-8d,%-4d,%-2d M\n", timems*1000/sum, HostMem.num_terminated_neutrons, 1, num_src/1000000);
  fclose(fp);
}

#if defined(__XS_GPU)
#include "CComplex.h"
void eval_xs(multipole mp_para,unsigned int* iS_h,unsigned int* iS_d, unsigned numIso, CMPTYPE E, CMPTYPE sqrtKT, 
                               CMPTYPE *sigTs_h, CMPTYPE *sigAs_h, CMPTYPE *sigFs_h,
                               CMPTYPE *sigTs_d, CMPTYPE *sigAs_d, CMPTYPE *sigFs_d){
  gpuErrchk(cudaMemcpy(iS_d, iS_h, sizeof(unsigned int), cudaMemcpyHostToDevice));  
  device_xs_eval<<<1,numIso>>>(mp_para,iS_d,E,sqrtKT,sigTs_d, sigAs_d, sigFs_d);
  gpuErrchk(cudaMemcpy(sigTs_h, sigTs_d, sizeof(CMPTYPE), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(sigAs_h, sigAs_d, sizeof(CMPTYPE), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(sigFs_h, sigFs_d, sizeof(CMPTYPE), cudaMemcpyDeviceToHost));
}

void printdevice(){
  cudaDeviceProp prop; 
  int count;
  cudaGetDeviceCount(&count);
  printf("num of devices=%d\n",count);
  for (int i=0; i<count; i++){
    cudaGetDeviceProperties( &prop, i );
    printf( "   --- General Information for device %d ---\n", i );
    printf( "Name:  %s\n", prop.name );
    printf( "Compute capability:  %d.%d\n", prop.major, prop.minor );
    printf( "Clock rate:  %d\n", prop.clockRate );
    printf( "Device copy overlap:  " );
    if (prop.deviceOverlap)
      printf( "Enabled\n" );
    else
      printf( "Disabled\n");
    printf( "Kernel execution timeout :  " );
    if (prop.kernelExecTimeoutEnabled)
      printf( "Enabled\n" );
    else
      printf( "Disabled\n" );
    
    printf( "   --- Memory Information for device %d ---\n", i );
    printf( "Total global mem:  %ld\n", prop.totalGlobalMem );
    printf( "Total constant Mem:  %ld\n", prop.totalConstMem );
    printf( "Max mem pitch:  %ld\n", prop.memPitch );
    printf( "Texture Alignment:  %ld\n", prop.textureAlignment );
    
    printf( "   --- MP Information for device %d ---\n", i );
    printf( "Multiprocessor count:  %d\n",
	    prop.multiProcessorCount );
    printf( "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
    printf( "Registers per mp:  %d\n", prop.regsPerBlock );
    printf( "Threads in warp:  %d\n", prop.warpSize );
    printf( "Max threads per block:  %d\n",
	    prop.maxThreadsPerBlock );
    printf( "Max thread dimensions:  (%d, %d, %d)\n",
	    prop.maxThreadsDim[0], prop.maxThreadsDim[1],
	    prop.maxThreadsDim[2] );
    printf( "Max grid dimensions:  (%d, %d, %d)\n",
	    prop.maxGridSize[0], prop.maxGridSize[1],
	    prop.maxGridSize[2] );
    printf( "\n" );
  }
}
#endif

#if defined(__W__GPU)
void eval_w(CPUComplex<CMPTYPE>* z_h, void** z_d, CPUComplex<CMPTYPE>* w_h, void** w_d,unsigned window){
  gpuErrchk(cudaMemcpy((CComplex<CMPTYPE>*)(*z_d),z_h,sizeof(CMPTYPE)*2*window,cudaMemcpyHostToDevice));
  device_w_eval<<<1,window>>>((CComplex<CMPTYPE>*)(*z_d),(CComplex<CMPTYPE>*)(*w_d)); 
  gpuErrchk(cudaMemcpy(w_h,(CComplex<CMPTYPE>*)(*w_d),sizeof(CMPTYPE)*2*window,cudaMemcpyDeviceToHost));
}
#endif

