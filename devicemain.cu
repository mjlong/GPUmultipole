#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"
#include "simulation.h"

/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/


void printdevice();

void anyvalue(struct multipoledata data, int setgridx, int setblockx){
  unsigned gridx, blockx, gridsize;
  unsigned ints=0, doubles=0, sharedmem;
  float timems = 0.0;
  unsigned *cnt;
  double *hostarray, *devicearray;
  struct neutronInfo Info;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  // printdevice();
  gridx = setgridx;
  blockx = setblockx;
  dim3 dimBlock(gridx, 1);
  dim3 dimGrid(blockx, 1, 1);
  gridsize = gridx*blockx;
  cudaMalloc((void**)&devicearray, 4*gridsize*sizeof(double));
  cudaMalloc((void**)&(Info.rndState), gridsize*sizeof(curandState));
  cudaMalloc((void**)&(Info.energy), gridsize*sizeof(double));
  cudaMalloc((void**)&(Info.tally), gridx*sizeof(unsigned));
  hostarray = (double*)malloc(4*gridsize*sizeof(double));
  cnt      = (unsigned*)malloc(gridx*sizeof(unsigned));

  multipole U238(data); //host multipoledata to device
  initialize<<<dimBlock, dimGrid>>>(Info, 2000.0);//1.95093e4);
  //  cudaDeviceSynchronize();
  /*
    Note: shared memory size is in unit of Bybe
    And the address can be referred in form of p = pshared + offset
  */
  ints = blockx;
  sharedmem = doubles*sizeof(double)+ints*sizeof(int);
  cudaEventRecord(start, 0);
  history<<<dimBlock, dimGrid, sharedmem>>>(U238, devicearray, Info);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&timems, start, stop);

  printf("time elapsed:%3.1f ms\n", timems);
 
  cudaMemcpy(hostarray, devicearray, 4*gridsize*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(cnt, Info.tally, gridx*sizeof(unsigned), cudaMemcpyDeviceToHost);

  for(int i=0;i<gridsize;i++){
    printf("%.15e %.15e %.15e %.15e\n",
	   hostarray[4*i],
	   hostarray[4*i+1],
	   hostarray[4*i+2],
	   hostarray[4*i+3]);
  }

  unsigned sum = 0;
  for (int i=0;i<gridx;i++){
    //printf("%4d\n",cnt[i]);
    sum += cnt[i];
  }
  printf("time elapsed:%g mus\n", timems*1000/sum);

  FILE *fp=NULL;
  fp = fopen("timelog","a+");
  fprintf(fp,"%3d,%3d,%g    \n", gridx, blockx, timems*1000/sum);
  fclose(fp);
  //cudaEventRecord(stop, 0);
  //cudaEventSynchronize(stop);
  //cudaEventElapsedTime(&timems, start, stop);

  //cudaEventDestroy(start);
  //cudaEventDestroy(stop);

  cudaFree(devicearray);
  cudaFree(Info.energy);
  cudaFree(Info.tally);
  cudaFree(Info.rndState);

  free(hostarray);
  free(cnt);
  return;
}


void printdevice(){
  cudaDeviceProp prop; 
  int count;
  cudaGetDeviceCount(&count);
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
