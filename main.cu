#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"
#include <cuda.h>
#include <curand_kernel.h>

/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/

__global__ void history(multipole, curandState *rndState, double *, double *);
void printdevice();

void anyvalue(struct multipoledata data, int *value, double *d1, double *d2){
  curandState *rndState;
  unsigned gridx, gridy, blockx, blocky, blockz, blocknum, gridsize;
  unsigned ints=0, floats=0, doubles=0, sharedmem;
  double *hostarray, *devicearray, *tally, *dev_tally;
  printdevice();
  gridx = 4;
  gridy = 4;
  blockx = 32;
  blocky = 1;
  blockz = 1;
  dim3 dimBlock(gridx, gridy);
  dim3 dimGrid(blockx, blocky, blockz);
  blocknum = gridx*gridy; 
  gridsize = gridx*gridy*blockx*blocky*blockz;
  cudaMalloc((void**)&rndState, gridsize*sizeof(curandState));
  cudaMalloc((void**)&devicearray, 7*gridsize*sizeof(double));
  cudaMalloc((void**)&dev_tally, blocknum*sizeof(double));
  hostarray = (double*)malloc(7*gridsize*sizeof(double));
  tally     = (double*)malloc(blocknum*sizeof(double));
  multipole U238(data); //host multipoledata to device
  
  /*
    Note: shared memory size is in unit of Bybe
    And the address can be referred in form of p = pshared + offset
  */
  doubles = blockx*blocky*blockz;
  sharedmem = doubles*sizeof(double)+floats*sizeof(float)+ints*sizeof(int);
  history<<<dimBlock, dimGrid, sharedmem>>>(U238, rndState, devicearray, dev_tally);
  
  cudaMemcpy(hostarray, devicearray, 7*gridsize*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(tally, dev_tally, blocknum*sizeof(double), cudaMemcpyDeviceToHost);

  for(int i=0;i<gridsize;i++){
    printf("%8.4f %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n",
	   hostarray[7*i],
	   hostarray[7*i+1],
	   hostarray[7*i+2],
	   hostarray[7*i+3],
	   hostarray[7*i+4],
	   hostarray[7*i+5],
	   hostarray[7*i+6]);
  }
  for (int i=0;i<blocknum;i++)
    printf("%2.1f\n",tally[i]);


  return;
}


__global__ void history(multipole U238, curandState *rndState, double *devicearray, double *dev_tally){
  int blocksize = blockDim.x * blockDim.y * blockDim.z;
  
  //size of shared[] is given as 3rd parameter while launching the kernel
  extern __shared__ double shared[];
  
  double *tally = &shared[0];

  int i;

  int idl = 
    (blockDim.x*blockDim.y)*threadIdx.z+
    blockDim.x*threadIdx.y+
    threadIdx.x;

  int idb = 
    blockIdx.y*gridDim.x+blockIdx.x;

  //TODO:this is one scheme to match threads to 1D array, 
  //try others when real simulation structure becomes clear
  int id = 
    (blockDim.x*blockDim.y*blockDim.z)*
    (blockIdx.y*gridDim.x+blockIdx.x)+
    (blockDim.x*blockDim.y)*threadIdx.z+
    blockDim.x*threadIdx.y+
    threadIdx.x;



  
  /* Each thread gets same seed, a different sequence number, no offset */
  curand_init(1234, id, 0, &rndState[id]);
  /* Copy state to local memory for efficiency */ 
  curandState localState = rndState[id];
  bool live=true;
  double energy = 1000.0;
  double rnd;
  double sigT, sigA, sigF, sibT, sibA, sibF;
  while(live){
    rnd = curand_uniform(&localState);
    energy = energy * rnd;
    //for test
    energy = (1+id)*1.63;
    U238.xs_eval_fast(energy, sigT, sigA, sigF);
    U238.xs_eval_fast(energy, sqrt(900*KB), sibT, sibA, sibF);

    live = false;//(energy>50.0);
  }
  
  devicearray[7*id]=energy;
  devicearray[7*id+1]=sibT;
  devicearray[7*id+2]=sigT;
  devicearray[7*id+3]=sibA;
  devicearray[7*id+4]=sigA;
  devicearray[7*id+5]=sibF;
  devicearray[7*id+6]=sigF;

  /* Copy state back to global memory */ 
  rndState[id] = localState; 

  /*reduce tally*/

  tally[idl] = (double)(rnd<0.5);
  __syncthreads();
  i = blocksize>>2;
  while(0!=i){
    if(idl<i)
      tally[idl] += tally[idl+i];
    __syncthreads();
    i=i>>1;
  }
  if(0==idl)
    dev_tally[idb] = tally[0]/blocksize;

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
