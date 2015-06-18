#include "simulation.h"
#include "manmemory.h"

#include "devicebridge.h"
/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/

void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int banksize,int ubat){
  srand(100);
  int i=0;
  for(i=0;i<ubat;i++){
  //  printf("init... %d:%d/%d\n",i*gridx*blockx,(i+1)*gridx*blockx,banksize);
    initialize<<<gridx, blockx>>>(DeviceMem,width,banksize,i*gridx*blockx);
  }
  //gpuErrchk(cudaDeviceSynchronize());  
#if defined(__3D)
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridx*blockx*ubat,DeviceMem.nInfo.pos_x,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridx*blockx*ubat,DeviceMem.nInfo.pos_y,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridx*blockx*ubat,DeviceMem.nInfo.pos_z,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
#endif
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
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize){
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
#if defined(__3D)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int num_srcp, int num_src, int csize, int ibat, int nbat){
  float* x2 = (float*)malloc(sizeof(float)*num_src*2);
  float* y2 = (float*)malloc(sizeof(float)*num_src*2);
  float* z2 = (float*)malloc(sizeof(float)*num_src*2);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*num_src);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live ,sizeof(int)*num_src,   cudaMemcpyDeviceToHost));  
  int live;  unsigned j=0;int k=0; int igen; int ic=0; int avastart = 0;
  /*
  for(int i=0;i<num_src;i++){
    printf("%d ",HostMem.nInfo.live[i]);
    if(0==i%100) printf("\n");
  }
  printf("\n");
  */
  for(int i=0;i<num_src;i++){
    ic = avastart;
    live = HostMem.nInfo.live[i];
    if(live>10){
      igen = ceil(rand()*1.0/RAND_MAX*HostMem.wdspp[7])+1+ibat; //delayed generations>=2
      if(igen<nbat){
        HostMem.newly_delayed[igen]+=1;
        while( (ic<csize)&&( HostMem.nInfo.d_igen[ic]>=ibat) ){
	  ic++;
	}//end while ic
        avastart = ic;
        if(csize>ic){
          HostMem.nInfo.d_igen[ic]  = igen;
          HostMem.nInfo.d_pos_x[ic] = HostMem.nInfo.pos_x[i];
          HostMem.nInfo.d_pos_y[ic] = HostMem.nInfo.pos_y[i];
          HostMem.nInfo.d_pos_z[ic] = HostMem.nInfo.pos_z[i];
          HostMem.nInfo.d_nu[ic]    = live/10;
        }
        else
          printf("[Warning]: live=%d,ic=%d(>csize=%d),insufficient memory for newly delayed neutrons\n",live,ic,csize);
      }//end if fissioned generation number is within range
    }//end if live = 20 or 30
    else{
      for(k=0;k<live;k++){//live=2 or 3
	if(j>(num_src*2)) {printf("live=%d,j=%d,i=%d/%d,overflow\n",live,j,i,num_src);exit(-1);}
	x2[j]=HostMem.nInfo.pos_x[i];
	y2[j]=HostMem.nInfo.pos_y[i];
	z2[j]=HostMem.nInfo.pos_z[i];
	j++;
      }
    }//end if live=2 or 3
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+num_srcp,x2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+num_srcp,y2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+num_srcp,z2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  free(x2);  free(y2);  free(z2);
  return j;
}
//In all setbank flavors, "gridsize" actually denotes gridsize*ubat(*factor)
unsigned setbank_prompt(MemStruct DeviceMem, MemStruct HostMem, int num_src){
  float* x2 = (float*)malloc(sizeof(float)*num_src*2);
  float* y2 = (float*)malloc(sizeof(float)*num_src*2);
  float* z2 = (float*)malloc(sizeof(float)*num_src*2);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*num_src, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*num_src);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live ,sizeof(int)*num_src,   cudaMemcpyDeviceToHost));  
  int live;  unsigned j=0;int k=0; 

  for(int i=0;i<num_src;i++){
    live = HostMem.nInfo.live[i];
    for(k=0;k<live;k++){//live=2 or 3
      if(j>(num_src*2)) {printf("live=%d,j=%d,i=%d/%d,overflow\n",live,j,i,num_src);exit(-1);}
      x2[j]=HostMem.nInfo.pos_x[i];
      y2[j]=HostMem.nInfo.pos_y[i];
      z2[j]=HostMem.nInfo.pos_z[i];
      j++;
    }
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+num_src,x2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+num_src,y2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+num_src,z2,sizeof(float)*num_src*2, cudaMemcpyHostToDevice));  
  free(x2);  free(y2);  free(z2);
  return j;
}


int add_delayed(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize, int csize, int ibat, int nbat, int banksize){
  //============================================================================
  //=============== Add new delayed neutrons ===================================
  int ic,igen,livi,j,unlivestart;
  ic=0; igen=0; //igen now used as counter for newly delayed neutrons already added
  unlivestart = 0; j=0;
  //printf("[%3d]:adding %d...:",ibat,HostMem.newly_delayed[ibat]);
  while( (igen<(HostMem.newly_delayed[ibat]))&&(ic<csize) ){
    //printf("(ic=%d,desgen=%d,live=%d)\n",ic,HostMem.nInfo.d_igen[ic],HostMem.nInfo.d_nu[ic]);
    if(ibat==HostMem.nInfo.d_igen[ic]){
      igen++;
      livi = HostMem.nInfo.d_nu[ic];
      while(0<livi){
        j = unlivestart;
        if(j>=gridsize){
          printf("error bank overflow\n");
          exit(-1);
        }
        unlivestart = j+1; //update unlive start
        HostMem.nInfo.pos_x[j] = HostMem.nInfo.d_pos_x[ic];
        HostMem.nInfo.pos_y[j] = HostMem.nInfo.d_pos_y[ic];
        HostMem.nInfo.pos_z[j] = HostMem.nInfo.d_pos_z[ic];
        livi-=1;
      }//end assigning
    }  //end if there exist neutron delayed at this generation
    ic++;
  }    //end search in delayed bank
  //============================================================================
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridsize+banksize, HostMem.nInfo.pos_x, sizeof(float)*j, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridsize+banksize, HostMem.nInfo.pos_y, sizeof(float)*j, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridsize+banksize, HostMem.nInfo.pos_z, sizeof(float)*j, cudaMemcpyHostToDevice));
  return j;
}


#endif


int count_pop(int *live, int gridsize){
  int sum = 0;
  for(int i=0;i<gridsize;i++)
    sum += (0!=live[i]);
  return sum;
}

#if defined(__3D)
void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned separator,unsigned banksize, unsigned isTally){
  int i=0;
  for(i=0;i<ubat;i++){//num_src is important as loop index, but useless in history<<<>>>
    //printf("i=%d/%d\n",i,num_src/(gridx*blockx));
    history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, separator,i*gridx*blockx,banksize,isTally);
  }
}

void check(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat){
  int i=0;
  printf("start of check\n");
  for(i=0;i<ubat;i++){
    preview_live<<<gridx, blockx>>>(DeviceMem, i*gridx*blockx);
  }
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
