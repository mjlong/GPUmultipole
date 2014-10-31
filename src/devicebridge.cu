#include "CPUComplex.h"
#include "CComplex.h"
#include "simulation.h"
#include "manmemory.h"

#include "devicebridge.h"
/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/

void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem){
  initialize<<<gridx, blockx>>>(DeviceMem, STARTENE);//1.95093e4);
}

void start_neutrons(unsigned gridx, unsigned blockx, unsigned numIsos, multipole mp_data, CMPTYPE* devicearray, MemStruct DeviceMem, unsigned num_src){
#if defined(__TRACK)
    history<<<gridx, blockx, blockx*sizeof(unsigned)>>>(numIsos, mp_data, devicearray, DeviceMem, num_src);
#else
    history<<<gridx, blockx, blockx*sizeof(unsigned)>>>(numIsos, mp_data, DeviceMem, num_src);
#endif
} 

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src){
  unsigned active;
  statistics<<<1, gridx, gridx*sizeof(unsigned)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_terminated_neutrons);
  gpuErrchk(cudaMemcpy(HostMem.num_terminated_neutrons,DeviceMem.num_terminated_neutrons,sizeof(unsigned int), cudaMemcpyDeviceToHost));
  active = HostMem.num_terminated_neutrons[0] + gridx*blockx < num_src;  
  return active;
}

void remain_neutrons(unsigned gridx, unsigned blockx, unsigned numIsos, multipole mp_data, CMPTYPE* devicearray, MemStruct DeviceMem){
  remaining<<<gridx, blockx>>>(numIsos, mp_data, devicearray, DeviceMem);
}

void print_results(unsigned gridx, unsigned blockx, unsigned num_src, MemStruct DeviceMem, MemStruct HostMem, CMPTYPE* hostarray, CMPTYPE* devicearray, unsigned* blockcnt,unsigned* cnt, float timems){
  gpuErrchk(cudaMemcpy(hostarray, devicearray, 4*gridx*blockx*sizeof(CMPTYPE), cudaMemcpyDeviceToHost));
  
  statistics<<<gridx, blockx, blockx*sizeof(int)>>>(DeviceMem.tally.cnt, blockcnt);
  gpuErrchk(cudaMemcpy(cnt, blockcnt, gridx*sizeof(unsigned), cudaMemcpyDeviceToHost));

/*print energy & XS (energies for __TRACK)*/
#if !defined(__PLOT)
  for(int i=0;i<gridx*blockx;i++){
    printf(" %.15e %.15e %.15e %.15e",
	   hostarray[4*i],
	   hostarray[4*i+1],
	   hostarray[4*i+2],
	   hostarray[4*i+3]);
    if(hostarray[4*i]<0)
      printf("error-:%d \n",i);
    else{
      if(hostarray[4*i]>20000.0)
	printf("error+:%d \n",i);
      else
	printf("\n");
    }
  }
#endif

/*print collision cnt and time*/
#if !defined(__PROCESS) && !defined(__TRACK) && !defined(__PLOT)
  unsigned sum = 0;
  for (int i=0;i<gridx;i++){
    printf("%4d\n",cnt[i]);
    sum += cnt[i];
  }
  printf("time elapsed:%g mus\n", timems*1000/sum);

  FILE *fp=NULL;
  fp = fopen("timelog","a+");
  gpuErrchk(cudaMemcpy(HostMem.num_terminated_neutrons, 
		       DeviceMem.num_terminated_neutrons, 
		       sizeof(unsigned int), 
		       cudaMemcpyDeviceToHost));
  fprintf(fp,"%-4d,%-4d,%-.6f,%-8d,%-4d,%-2d M\n", gridx, blockx,timems*1000/sum, *HostMem.num_terminated_neutrons, 1, num_src/1000000);
  fclose(fp);
#endif
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
