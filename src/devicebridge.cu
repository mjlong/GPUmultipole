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
  initialize<<<gridx, blockx>>>(DeviceMem);
}

void start_neutrons(unsigned gridx, unsigned blockx, material mat, multipole mp_data, MemStruct DeviceMem, unsigned num_src){
    history<<<gridx, blockx, blockx*sizeof(unsigned)>>>(mat, mp_data, DeviceMem, num_src);
} 

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src){
//count terminated neutrons 
  unsigned active;
  reduce_sum_plus<<<1, gridx, gridx*sizeof(unsigned)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_terminated_neutrons);
  gpuErrchk(cudaMemcpy(HostMem.num_terminated_neutrons,DeviceMem.num_terminated_neutrons,sizeof(unsigned int), cudaMemcpyDeviceToHost));
  active = HostMem.num_terminated_neutrons[0] + gridx*blockx < num_src;  
  return active;
}

unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem){
//count neutrons still marked "live"
  unsigned active;
  reduce_sum_plus<<<1, gridx, gridx*sizeof(unsigned)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_terminated_neutrons);
  reduce_sum_equal<<<gridx,blockx,blockx*sizeof(unsigned)>>>(DeviceMem.nInfo.live, DeviceMem.block_terminated_neutrons);
  reduce_sum_equal<<<1,gridx, gridx*sizeof(unsigned)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_live_neutrons);
  gpuErrchk(cudaMemcpy(&active, DeviceMem.num_live_neutrons, sizeof(unsigned), cudaMemcpyDeviceToHost));  
  return active;
}

void sort_prepare(unsigned gridx, unsigned blockx,MemStruct DeviceMem, material mat){
  update_sort_key<<<gridx, blockx>>>(DeviceMem, mat);
}

void transport_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem, material mat, unsigned renew){
  transport<<<gridx, blockx>>>(DeviceMem, mat,renew);
}

void print_results(unsigned gridx, unsigned blockx, unsigned num_src, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem, unsigned* blockcnt,unsigned* cnt, float timems){
  
  unsigned *d_cnt, *h_cnt;
  gpuErrchk(cudaMalloc((void**)&d_cnt, num_bin*sizeof(unsigned)));
  h_cnt = (unsigned*)malloc(num_bin*sizeof(unsigned));
  for(int i=0;i<num_bin;i++){
    reduce_sum_plus<<<gridx, blockx, blockx*sizeof(int)>>>(DeviceMem.tally.cnt+i*gridx*blockx, blockcnt+i*gridx);
  }
  for(int i=0;i<num_bin;i++){
    reduce_sum_equal<<<1, gridx, gridx*sizeof(int)>>>(blockcnt+i*gridx, d_cnt+i);
  }
  gpuErrchk(cudaMemcpy(h_cnt,d_cnt,sizeof(unsigned)*num_bin, cudaMemcpyDeviceToHost));

  gpuErrchk(cudaMemcpy(cnt, blockcnt, num_bin*gridx*sizeof(unsigned), cudaMemcpyDeviceToHost));

/*print collision cnt and time*/
  unsigned sum = 0;
  for (int i=0;i<gridx;i++){
    for(int j=0;j<num_bin;j++){ 
      printf("%4d ",cnt[i+j*gridx]);
      sum += cnt[i+j*gridx];
    }
    printf("\n");
  }
  for(int j=0;j<num_bin;j++){ 
    printf("%4d ",h_cnt[j]);
  }
  printf("\n");
  printf("time elapsed:%g mus\n", timems*1000/sum);
  
  free(h_cnt);
  gpuErrchk(cudaFree(d_cnt));
  FILE *fp=NULL;
  fp = fopen("timelog","a+");
  gpuErrchk(cudaMemcpy(HostMem.num_terminated_neutrons, 
		       DeviceMem.num_terminated_neutrons, 
		       sizeof(unsigned int), 
		       cudaMemcpyDeviceToHost));
  fprintf(fp,"%-4d,%-4d,%-.6f,%-8d,%-4d,%-2d M\n", gridx, blockx,timems*1000/sum, *HostMem.num_terminated_neutrons, 1, num_src/1000000);
  fclose(fp);
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
