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


#if defined (__QUICKW)
#include "QuickW.h"
#endif

#if defined (__FOURIERW)
#include "fourierw.h"
__constant__ CMPTYPE a[M+1];
__constant__ CMPTYPE b[M+1];
#endif

#if defined (__QUICKWC) || defined(__INTERPEXP)
__constant__ CMPTYPE2 constwtable[LENGTH*LENGTH];
#endif


void printdevice();

void anyvalue(struct multipoledata data, unsigned setgridx, unsigned setblockx, unsigned num_src, unsigned devstep){
  unsigned gridx, blockx, gridsize;
  unsigned ints=0, sharedmem;
  float timems = 0.0;
  unsigned *cnt, *blockcnt;
  unsigned int active,i;
  CMPTYPE *hostarray, *devicearray;
  MemStruct HostMem, DeviceMem;
  cudaEvent_t start, stop;
  //printdevice();
  gpuErrchk(cudaEventCreate(&start));
  gpuErrchk(cudaEventCreate(&stop));
  gridx = setgridx;
  blockx = setblockx;
  dim3 dimBlock(gridx, 1);
  dim3 dimGrid(blockx, 1, 1);
  gridsize = gridx*blockx;
  gpuErrchk(cudaMalloc((void**)&devicearray, 4*gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMemset(devicearray, 0, 4*gridsize*sizeof(CMPTYPE)));
  gpuErrchk(cudaMalloc((void**)&(DeviceMem.nInfo), gridsize*sizeof(NeutronInfoStruct)));
  gpuErrchk(cudaMalloc((void**)&(DeviceMem.thread_active), gridsize*sizeof(unsigned int)));
  HostMem.thread_active = (unsigned int *)malloc(gridsize*sizeof(unsigned int));
  gpuErrchk(cudaMalloc((void**)&(DeviceMem.num_terminated_neutrons), sizeof(unsigned int)));
  HostMem.num_terminated_neutrons = (unsigned int *)malloc(sizeof(unsigned int));
  HostMem.num_terminated_neutrons[0] = 0u;
  gpuErrchk(cudaMemcpy(DeviceMem.num_terminated_neutrons, HostMem.num_terminated_neutrons, sizeof(unsigned int), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc((void**)&(DeviceMem.tally), gridsize*sizeof(TallyStruct)));
  gpuErrchk(cudaMalloc((void**)&(blockcnt), gridx*sizeof(unsigned int)));
  hostarray = (CMPTYPE*)malloc(4*gridsize*sizeof(CMPTYPE));
  cnt      = (unsigned*)malloc(gridx*sizeof(unsigned));

// construct coefficients a[n] for fourier expansion w
#if defined(__FOURIERW)
  CMPTYPE *da;
  CMPTYPE *db;
  gpuErrchk(cudaMalloc((void**)&da, (M+1)*sizeof(CMPTYPE))); 
  gpuErrchk(cudaMalloc((void**)&db, (M+1)*sizeof(CMPTYPE))); 
  fill_a<<<1,M+1>>>(da,db); 
  cudaMemcpyToSymbol(a, da, M*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
  cudaMemcpyToSymbol(b, db, M*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
#endif

// fill w function table for quickw
#if defined(__QUICKW)
  CComplex<CMPTYPE> *wtable;
  gpuErrchk(cudaMalloc((void**)&wtable, LENGTH*LENGTH * 2 * sizeof(CMPTYPE)));
  fill_w_tabulated<<<LENGTH,LENGTH>>>(wtable);
#if defined(__QUICKWC)
  cudaMemcpyToSymbol(constwtable, wtable, LENGTH*LENGTH*2*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
#endif
#if defined(__QUICKWT)
  bindwtable(wtable);
#endif
#endif

#if defined(__QUICKWG)
  multipole U238(data, wtable);
#else
  multipole U238(data);
#endif 

// fill exp(z) table for fourierw
#if defined(__INTERPEXP)
  CComplex<CMPTYPE> *exptable;
  gpuErrchk(cudaMalloc((void**)&exptable, LENGTH*LENGTH * 2 * sizeof(CMPTYPE)));
  fill_exp_table<<<LENGTH,LENGTH>>>(exptable);
  cudaMemcpyToSymbol(constwtable, exptable, LENGTH*LENGTH*2*sizeof(CMPTYPE), 0, cudaMemcpyDeviceToDevice);
#endif
  initialize<<<dimBlock, dimGrid>>>(DeviceMem, STARTENE);//1.95093e4);
  //  cudaDeviceSynchronize();
  /*
    Note: shared memory size is in unit of Bybe
    And the address can be referred in form of p = pshared + offset
  */
  gpuErrchk(cudaEventRecord(start, 0));

#if defined(__PROCESS) //|| defined(__TRACK)
  active = 0u;
#else
  active = 1u;
#endif

  while (active){
#if defined(__TRACK)
    history<<<dimBlock, dimGrid>>>(U238, devicearray, DeviceMem, num_src, devstep);
#else
    history<<<dimBlock, dimGrid>>>(U238, DeviceMem, num_src, devstep);
#endif
    gpuErrchk(cudaMemcpy(HostMem.thread_active, DeviceMem.thread_active, gridsize*sizeof(unsigned int), cudaMemcpyDeviceToHost));
    active = 0u;
    for (i = 0; i < blockx; i++){
      active += HostMem.thread_active[i];
    }
  }

  remaining<<<dimBlock, dimGrid>>>(U238, devicearray, DeviceMem);

  gpuErrchk(cudaEventRecord(stop, 0));
  gpuErrchk(cudaEventSynchronize(stop));
  gpuErrchk(cudaEventElapsedTime(&timems, start, stop));

  printf("time elapsed:%3.1f ms\n", timems);
 
  gpuErrchk(cudaMemcpy(hostarray, devicearray, 4*gridsize*sizeof(CMPTYPE), cudaMemcpyDeviceToHost));

  
  ints = blockx;
  sharedmem = ints*sizeof(int);
  statistics<<<dimBlock, dimGrid, sharedmem>>>(DeviceMem.tally, blockcnt);
  gpuErrchk(cudaMemcpy(cnt, blockcnt, gridx*sizeof(unsigned), cudaMemcpyDeviceToHost));

/*print energy & XS (energies for __TRACK)*/
#if !defined(__PLOT)
  for(int i=0;i<gridsize;i++){
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
  fprintf(fp,"%-4d,%-4d,%-.6f,%-8d,%-4d,%-2d M\n", gridx, blockx,timems*1000/sum, *HostMem.num_terminated_neutrons, devstep, num_src/1000000);
  fclose(fp);
#endif
  //cudaEventRecord(stop, 0);
  //cudaEventSynchronize(stop);
  //cudaEventElapsedTime(&timems, start, stop);

  gpuErrchk(cudaEventDestroy(start));
  gpuErrchk(cudaEventDestroy(stop));

  gpuErrchk(cudaFree(devicearray));
  gpuErrchk(cudaFree(DeviceMem.nInfo));
  gpuErrchk(cudaFree(DeviceMem.thread_active));
  gpuErrchk(cudaFree(DeviceMem.tally));
#if defined(__QUICKW)
  gpuErrchk(cudaFree(wtable));
#endif
#if defined(__FOURIERW)
  gpuErrchk(cudaFree(da));
  gpuErrchk(cudaFree(db));
#endif
  U238.release_pointer();

  free(hostarray);
  free(cnt);
  free(HostMem.thread_active);
  free(HostMem.num_terminated_neutrons);
  return;
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
