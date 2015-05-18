#include "simulation.h"
#include "manmemory.h"

#include "devicebridge.h"
/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/

void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int banksize){
  initialize<<<gridx, blockx>>>(DeviceMem,width,banksize);
}

#if defined(__SCATTERPLOT)
void copyinitial(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize){
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*gridsize, cudaMemcpyDeviceToHost)); 
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live,DeviceMem.nInfo.live,sizeof(int)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.energy,DeviceMem.nInfo.energy,sizeof(CMPTYPE)*gridsize, cudaMemcpyDeviceToHost));  
}
#endif

void resetcount(MemStruct DeviceMem){
  unsigned x=0;
  gpuErrchk(cudaMemcpy(DeviceMem.num_terminated_neutrons,&x,sizeof(unsigned), cudaMemcpyHostToDevice));  
}
#if defined(__1D)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize){
  float* y2 = (float*)malloc(sizeof(float)*gridsize);
  float* x2 = (float*)malloc(sizeof(float)*gridsize*3);
  gpuErrchk(cudaMemcpy(y2,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  float y;
  unsigned j=0;
  for(int i=0;i<gridsize;i++){
    y = y2[i];
    if(0!=y){
      if(y>0){
	//number=3;
	x2[j++]=y;
	x2[j++]=y;
	x2[j++]=y;
      }
      else{
	//number=2;
	x2[j++]=0-y;
	x2[j++]=0-y;
      }
    }
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x,x2,sizeof(float)*gridsize*3, cudaMemcpyHostToDevice));  
  free(x2);
  free(y2);
  return j;
}
#endif
#if defined(__3D)&&!defined(__TRAN)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize){
  float* x2 = (float*)malloc(sizeof(float)*gridsize*3);
  float* y2 = (float*)malloc(sizeof(float)*gridsize*3);
  float* z2 = (float*)malloc(sizeof(float)*gridsize*3);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live ,sizeof(int)*gridsize,   cudaMemcpyDeviceToHost));  
  int live;  unsigned j=0;int k=0;
  for(int i=0;i<gridsize;i++){
    live = HostMem.nInfo.live[i];
    for(k=0;k<live;k++){//live=2 or 3
      x2[j]=HostMem.nInfo.pos_x[i];
      y2[j]=HostMem.nInfo.pos_y[i];
      z2[j]=HostMem.nInfo.pos_z[i];
      j++;
    }
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x,x2,sizeof(float)*gridsize*3, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y,y2,sizeof(float)*gridsize*3, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z,z2,sizeof(float)*gridsize*3, cudaMemcpyHostToDevice));  
  free(x2);  free(y2);  free(z2);
  return j;
}
#endif

int flushbank(MemStruct DeviceMem, MemStruct HostMem,unsigned lastpop,float a,unsigned gridsize){
  //gridsize = num_src = factor*gridx*blockx
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*gridsize, cudaMemcpyDeviceToHost)); 
  gpuErrchk(cudaMemcpy(HostMem.nInfo.dir_polar,DeviceMem.nInfo.dir_polar,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.dir_azimu,DeviceMem.nInfo.dir_azimu,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.d_closest,DeviceMem.nInfo.d_closest,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live,DeviceMem.nInfo.live,sizeof(int)*gridsize, cudaMemcpyDeviceToHost));  
#if defined(__SCATTERPLOT)
  gpuErrchk(cudaMemcpy(HostMem.nInfo.energy,DeviceMem.nInfo.energy,sizeof(CMPTYPE)*gridsize, cudaMemcpyDeviceToHost));  
#endif
  //for(int i=0;i<gridsize;i++)
  //  printf("%2d ", HostMem.nInfo.live[i]);
  //printf("[l:%d]\n",lastpop);

  unsigned unlivestart=0;
  int i,j,ilp,inp,inp2,livi;
  ilp = 0; i=0; inp=0; inp2=0;
  while((ilp<lastpop)&&(i<gridsize)){// I assert ilp reaches lastpop no later than i reaches gridsize
    livi = HostMem.nInfo.live[i];
    //*allOld = (*allOld)&&((livi<=-3));

    ilp += (0!=livi)&&(-2!=livi)&&(-4!=livi); 
    inp += (1<=livi)*livi;
    inp2+= (1<=livi);
    //printf("i=%d, ilp=%d, lastpop=%d, live=%d\n",i,ilp,lastpop,livi);

    while(1<livi){
      //live=1 continue; 
      //live=0 didn't run; 
      //live=-1 terminated; 
      //live=-2 refreshed by host
      //live>1 fission to live neutrons
      j = unlivestart;
      while((1<=HostMem.nInfo.live[j])||(-1>HostMem.nInfo.live[j])){//live==-1 or 0 can be refreshed
	j++;
      }
      //printf("      live[%d] is changed from %d to -2\n",j,HostMem.nInfo.live[j]);
      if(j>=gridsize){
	printf("error bank overflow\n");
	exit(-1);
      }
      unlivestart = j+1; //update unlive start
      HostMem.nInfo.pos_x[j] = HostMem.nInfo.pos_x[i];
      HostMem.nInfo.pos_y[j] = HostMem.nInfo.pos_y[i];
      HostMem.nInfo.pos_z[j] = HostMem.nInfo.pos_z[i];
      HostMem.nInfo.dir_polar[j] = HostMem.nInfo.dir_polar[i];
      HostMem.nInfo.dir_azimu[j] = HostMem.nInfo.dir_azimu[i];
      HostMem.nInfo.d_closest[j] = HostMem.nInfo.d_closest[i];
#if defined(__SCATTERPLOT)
      HostMem.nInfo.energy[j] = HostMem.nInfo.energy[i];
#endif
      lastpop-=((j>i)&&(-1==HostMem.nInfo.live[j]));//if later j is reflushed, I don't want it to be counted in ilp
      HostMem.nInfo.live[j] = -2;

      inp2 += 1; //second next generation population counter
      livi-=1; HostMem.nInfo.live[i]-=1;
    }
    i++;
  }
  //for(int i=0;i<gridsize;i++)
  //  printf("%2d ", HostMem.nInfo.live[i]);
  //printf("[n:%d?=%d]\n",inp,inp2);

  if(0==inp){
    for(int i=0;i<gridsize;i++){
      livi = HostMem.nInfo.live[i];
      HostMem.nInfo.live[i] = (-3==livi)||(-4==livi);
      //printf("%2d ", HostMem.nInfo.live[i]);
    }
    //printf("\n");
  }
  //printf("\n");
  //If a threads has live=-1 or 0 but has not been refreshed here, it must be treated with care at first of history<<<>>>
  //all possible live are: -1 terminated and not refreshed; 
  //                        0 didn't run and not refreshed
  //                       -2 refreshed, need direction sample
  //                        1 continue
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x,    HostMem.nInfo.pos_x,    sizeof(float)*gridsize, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y,    HostMem.nInfo.pos_y,    sizeof(float)*gridsize, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z,    HostMem.nInfo.pos_z,    sizeof(float)*gridsize, cudaMemcpyHostToDevice)); 
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.dir_polar,HostMem.nInfo.dir_polar,sizeof(float)*gridsize, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.dir_azimu,HostMem.nInfo.dir_azimu,sizeof(float)*gridsize, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.d_closest,HostMem.nInfo.d_closest,sizeof(float)*gridsize, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.live,     HostMem.nInfo.live,     sizeof(int)  *gridsize, cudaMemcpyHostToDevice));  

  return inp;
}

int count_pop(int *live, int gridsize){
  int sum = 0;
  for(int i=0;i<gridsize;i++)
    sum += (1==live[i]);
  return sum;
}

#if defined(__3D)&&defined(__TRAN)
void transient_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize,float p2){
  int i;
  for(i=0;(i*gridx*blockx)<num_src;i++){
    //printf("i=%d/%d\n",i,num_src/(gridx*blockx));
    history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, num_src,active,i*gridx*blockx,p2);
  }
} 
#else
void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize){
  history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, num_src,active,banksize);
} 
#endif

//Original branches of start_neutron() for 1D,3D,ref,vac and steady, transient
//void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize){
//#if defined(__3D)&&!defined(__TRAN)
//  history_3d_ref<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, num_src,active,banksize);
//#endif
//} 
//
//
unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src){
//count terminated neutrons 
  unsigned active;
  reduce_sum_plus<<<1, gridx, gridx*sizeof(int)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_terminated_neutrons);
  gpuErrchk(cudaMemcpy(HostMem.num_terminated_neutrons,DeviceMem.num_terminated_neutrons,sizeof(int), cudaMemcpyDeviceToHost));
  active = HostMem.num_terminated_neutrons[0] + gridx*blockx < num_src;  
#if defined(__PRINTTRACK__)
  printf("[active]%d terminated\n",HostMem.num_terminated_neutrons[0]);
#endif
  return active;
}

unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem){
//count neutrons still marked "live"
  int active;
  reduce_sum_equal<<<gridx,blockx,blockx*sizeof(int)>>>(DeviceMem.nInfo.live, DeviceMem.block_terminated_neutrons);
  //I made a mistake to reuse block_terminated_neutrons here. 
  //However, as long as blockx<=gridx(size of block_terminated_neutrons), there would be no problem
  reduce_sum_equal<<<1,gridx, gridx*sizeof(int)>>>(DeviceMem.block_terminated_neutrons, DeviceMem.num_live_neutrons);
  gpuErrchk(cudaMemcpy(&active, DeviceMem.num_live_neutrons, sizeof(int), cudaMemcpyDeviceToHost));  
  return active;
}

void save_results(unsigned ibat, unsigned gridx, unsigned blockx, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem){
  for(int i=0;i<num_bin;i++){
    reduce_sum_equal<<<gridx, blockx, blockx*sizeof(CMPTYPE)>>>(
                   DeviceMem.tally.cnt+i*gridx*blockx, 
                   DeviceMem.block_spectrum+i*gridx);
  }
  for(int i=0;i<num_bin;i++){
    reduce_sum_equal<<<1, gridx, gridx*sizeof(CMPTYPE)>>>(
                   DeviceMem.block_spectrum+i*gridx, DeviceMem.batcnt+i);
  }
  //printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  //printf("%s\n", cudaGetErrorString(cudaThreadSynchronize()));
  gpuErrchk(cudaMemcpy(HostMem.batcnt,DeviceMem.batcnt,sizeof(CMPTYPE)*num_bin, cudaMemcpyDeviceToHost));

/*print collision cnt and time*/
/*
  unsigned sum=0;
  for(int j=0;j<num_bin;j++){ 
    sum+=h_cnt[j];
    printf("%6d ",h_cnt[j]);
  }
  printf("|||%u +++ %u\n",HostMem.num_terminated_neutrons[0],sum);
*/
}

void print_results(unsigned meshes, unsigned nbat, double *tally){
  int im,ib;
  for(ib=0;ib<nbat;ib++){
    for(im=0;im<meshes;im++){
      printf("%.5f ",tally[ib*meshes+im]);
    }
    printf("\n");
  }
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
