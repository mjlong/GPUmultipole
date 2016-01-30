#include "simulation.h"
#include "manmemory.h"

#include "devicebridge.h"
/*
  To compile host and device codes separately, 
  this "main" file works as interface 
  allocating device memory, transfering data and partitioning computation sources
*/
void initialize_neutrons_fix(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int ubat){
  srand (time(NULL));
  int i=0;
  for(i=0;i<ubat;i++){
    fixsrc_sample<<<gridx, blockx>>>(DeviceMem,width,i*gridx*blockx);
  }
  //gpuErrchk(cudaDeviceSynchronize());  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridx*blockx*ubat,DeviceMem.nInfo.pos_x,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridx*blockx*ubat,DeviceMem.nInfo.pos_y,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridx*blockx*ubat,DeviceMem.nInfo.pos_z,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
}

void initialize_neutrons_active_not_src(unsigned gridx, unsigned blockx,MemStruct DeviceMem, int num_seg, int seed){
  int i=0;
  for(i=0;i<num_seg;i++){
    initialize_without_src<<<gridx, blockx>>>(DeviceMem,i*gridx*blockx,seed);
  }
  //gpuErrchk(cudaDeviceSynchronize());  
}


void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int banksize,int ubat, int seed){
  int i=0;
  for(i=0;i<ubat;i++){
  //  printf("init... %d:%d/%d\n",i*gridx*blockx,(i+1)*gridx*blockx,banksize);
    initialize<<<gridx, blockx>>>(DeviceMem,width,banksize,i*gridx*blockx,seed);
  }
  //gpuErrchk(cudaDeviceSynchronize());  
#if defined(__3D)
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridx*blockx*ubat,DeviceMem.nInfo.pos_x,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridx*blockx*ubat,DeviceMem.nInfo.pos_y,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridx*blockx*ubat,DeviceMem.nInfo.pos_z,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
#endif
#if defined(__1D)
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridx*blockx*ubat,DeviceMem.nInfo.pos_x,sizeof(float)*gridx*blockx*ubat, cudaMemcpyDeviceToDevice));    
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
#if defined(__MTALLY)||(__FTALLY)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, int tnum_bins){
  float* y2 = (float*)malloc(sizeof(float)*gridsize);
  float* x2 = (float*)malloc(sizeof(float)*gridsize*2);
  int* sid1 = (int*)malloc(sizeof(int)*gridsize);
#if defined(__MTALLY)
  int* sid2 = (int*)malloc(sizeof(int)*gridsize*2);
#endif
  gpuErrchk(cudaMemcpy(y2,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(sid1,DeviceMem.nInfo.imat,sizeof(int )*gridsize, cudaMemcpyDeviceToHost));  
  int sid;
  float y; 
  unsigned j=0;
  for(int i=0;i<gridsize;i++){
    y = y2[i]; sid = sid1[i];
    HostMem.batcnt[sid]++;
    if(0!=y){
#if defined(__MTALLY)
      sid = sid/tnum_bins;
      if(y>0){sid2[j]=sid; x2[j++]=y;sid2[j]=sid;x2[j++]=y;sid2[j]=sid;x2[j++]=y;}
      else{sid2[j]=sid; x2[j++]=0-y;sid2[j]=sid; x2[j++]=0-y;}
#else
      if(y>0){x2[j++]=y;x2[j++]=y;x2[j++]=y;}
      else{x2[j++]=0-y;x2[j++]=0-y;}
#endif
    }
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridsize,x2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  free(x2);
  free(y2);
#if defined(__MTALLY)
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.imat+gridsize,sid2,sizeof(int)*gridsize*2, cudaMemcpyHostToDevice));  
  free(sid2);
#endif
  free(sid1);
  return j;
}
#else
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize){
  float* y2 = (float*)malloc(sizeof(float)*gridsize);
  float* x2 = (float*)malloc(sizeof(float)*gridsize*2);
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
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridsize,x2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  free(x2);
  free(y2);
  return j;
}
#endif//__MTALLY
#endif//__1D
#if defined(__3D)
//==============================================================================
//================setbank_converge() ===========================================
//1. Update fission sites in the phase of converging fission source
//2. The update follows traditional method, multiplicity is not treated
unsigned setbank_converge(MemStruct DeviceMem, MemStruct HostMem, int gridsize){
  float* x2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* y2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* z2 = (float*)malloc(sizeof(float)*gridsize*2);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live ,sizeof(int)*gridsize,   cudaMemcpyDeviceToHost));  
  int live;  unsigned j=0;int k=0;
  for(int i=0;i<gridsize;i++){
    live = HostMem.nInfo.live[i];
    //if(live<4){
    for(k=0;k<live;k++){//live=2 or 3
      if(j>(gridsize*2)) {printf("live=%d,j=%d,i=%d/%d,overflow\n",live,j,i,gridsize);exit(-1);}
      //else{
      x2[j]=HostMem.nInfo.pos_x[i];
      y2[j]=HostMem.nInfo.pos_y[i];
      z2[j]=HostMem.nInfo.pos_z[i];
      j++;
      //}
    }
    //}
  }
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridsize,x2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridsize,y2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridsize,z2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  free(x2);  free(y2);  free(z2);
  return j;
}
//=====================end function setbank_converge() =========================


//==============================================================================
//================setbank_prepbank() ===========================================
//1. Update fission sites in the phase of preparing delayed fission bank =======
//2. The update follows traditional method, multiplicity is not treated ========
unsigned setbank_prepbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, unsigned ibat){
  float* x2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* y2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* z2 = (float*)malloc(sizeof(float)*gridsize*2);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live ,sizeof(int)*gridsize,   cudaMemcpyDeviceToHost));  
  int live;  unsigned j=0;int k=0;
  unsigned cursor = HostMem.bank.cursor_end[0];
  for(int i=0;i<gridsize;i++){
    live = HostMem.nInfo.live[i]; 
    //==========If fissioned, fission site generates neutrons for the next generation ==============
    //======================= fission site also generates neutrons into the delay bank =============
    for(k=0;k<live;k++){//live=2 or 3
      if(j>(gridsize*2)) {printf("live=%d,j=%d,i=%d/%d,overflow\n",live,j,i,gridsize);exit(-1);}
      //else{
      x2[j]=HostMem.nInfo.pos_x[i];
      y2[j]=HostMem.nInfo.pos_y[i];
      z2[j]=HostMem.nInfo.pos_z[i];
      if(cursor<HostMem.bank.size[0]){
      HostMem.bank.x[cursor] = HostMem.nInfo.pos_x[i];
      HostMem.bank.y[cursor] = HostMem.nInfo.pos_y[i];
      HostMem.bank.z[cursor] = HostMem.nInfo.pos_z[i];
      HostMem.bank.available[cursor] = 1; 
      cursor++;
      }
      j++;
      //}
    }
  }
  HostMem.bank.cursor_end[0] = cursor;
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+gridsize,x2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+gridsize,y2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+gridsize,z2,sizeof(float)*gridsize*2, cudaMemcpyHostToDevice));  
  free(x2);    free(y2);    free(z2);
  return j;
}
//=====================end function setbank_prepbank() =========================


void setbank_active_out(MemStruct DeviceMem, MemStruct HostMem, unsigned banksize){
  float* x2 = (float*)malloc(sizeof(float)*(banksize));
  float* y2 = (float*)malloc(sizeof(float)*(banksize));
  float* z2 = (float*)malloc(sizeof(float)*(banksize));

  bank_pull(HostMem,x2,y2,z2,banksize); 

  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+banksize,x2,sizeof(float)*(banksize), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+banksize,y2,sizeof(float)*(banksize), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+banksize,z2,sizeof(float)*(banksize), cudaMemcpyHostToDevice));  

  free(x2);  free(y2);  free(z2);
}

void bank_pull(MemStruct HostMem, float *x2, float *y2, float* z2, unsigned banksize){
  int maxindex = HostMem.bank.cursor_end[0]; 
  int index=0;
  for(int i=0; i<banksize; i++){
    do{
      index = rand()%maxindex; 
    }while( 0==HostMem.bank.available[index] );

    x2[i] = HostMem.bank.x[index];
    y2[i] = HostMem.bank.y[index];
    z2[i] = HostMem.bank.z[index];
    //printf("-------------------------------pull %d to %d\n", index, i);
    HostMem.bank.available[index] = 0;
    HostMem.pull_list->EnQueue(index); //HostMem.pull_list->ViewQueue();
    if(HostMem.pull_list->IsFull()){
      printf("Warning: pull_list is full, there are bank positions that will never be used\n"); exit(-1);
    }

  }
}


void setbank_active_in(MemStruct DeviceMem, MemStruct HostMem, int gridsize, int banksize, int shift){
  int* fission_sites = (int*)malloc(sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(fission_sites,DeviceMem.nInfo.imat+shift,sizeof(int )*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live +shift,sizeof(int)*gridsize,   cudaMemcpyDeviceToHost));  
  int live;  int k=0; int fission_site;
  for(int i=0; i<gridsize; i++){
    live = HostMem.nInfo.live[i];
    fission_site = fission_sites[i];
    HostMem.batcnt[fission_site]+= (1*(0!=live)); 
    for(k=0;k<live;k++){
      bank_push(HostMem, i);
    }
  }
  gpuErrchk(cudaDeviceSynchronize());  
  free(fission_sites);
}

void bank_push(MemStruct HostMem, unsigned index){
  int j;
  if(!(HostMem.pull_list->IsEmpty()) ){
    j = HostMem.pull_list->DeQueue(); //HostMem.pull_list->ViewQueue();
    //printf("=======================pushing back to %d\n", j);
    HostMem.bank.x[j]=HostMem.nInfo.pos_x[index];
    HostMem.bank.y[j]=HostMem.nInfo.pos_y[index];
    HostMem.bank.z[j]=HostMem.nInfo.pos_z[index];
    HostMem.bank.available[j] = 1;
    return;
  }

  if(HostMem.bank.cursor_end[0]<HostMem.bank.size[0]){
    j = HostMem.bank.cursor_end[0];
    //printf("=======================PUSHing back to %d\n", j);
    HostMem.bank.x[j]=HostMem.nInfo.pos_x[index];
    HostMem.bank.y[j]=HostMem.nInfo.pos_y[index];
    HostMem.bank.z[j]=HostMem.nInfo.pos_z[index];
    HostMem.bank.available[j] = 1;
    HostMem.bank.cursor_end[0] = j+1;
    return;
  }
  
  printf("pusing %d, Queue is Empty and bank is full\n",index);
} 

#if defined(__FTALLY2)
//==============================================================================
//=============setbank() of __FTALLY2 samples to meet the \mu===================
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, int banksize, unsigned jstart, int shift){
  float* x2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* y2 = (float*)malloc(sizeof(float)*gridsize*2);
  float* z2 = (float*)malloc(sizeof(float)*gridsize*2);
  int* sid1 = (int*)malloc(sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(sid1,DeviceMem.nInfo.imat,sizeof(int )*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_x,DeviceMem.nInfo.pos_x+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_y,DeviceMem.nInfo.pos_y+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  gpuErrchk(cudaMemcpy(HostMem.nInfo.pos_z,DeviceMem.nInfo.pos_z+shift,sizeof(float)*gridsize, cudaMemcpyDeviceToHost));  
  memset(HostMem.nInfo.live,0,sizeof(int)*gridsize);
  gpuErrchk(cudaMemcpy(HostMem.nInfo.live, DeviceMem.nInfo.live +shift,sizeof(int)*gridsize,   cudaMemcpyDeviceToHost));  
  //jfor(int i=0;i<gridsize;i++){
  //j  printf("i=%d,live=%d,imat[i]=%d\n",i,HostMem.nInfo.live[i],sid1[i]);
  //j}
  int live;  unsigned j=jstart;int k=0; int sid;
  for(int i=0;i<gridsize;i++){
    live = HostMem.nInfo.live[i];
    sid = sid1[i];
    HostMem.batcnt[sid]+= (1*(0!=live));
    if(live>1){
    x2[j-jstart]=HostMem.nInfo.pos_x[i];
    y2[j-jstart]=HostMem.nInfo.pos_y[i];
    z2[j-jstart]=HostMem.nInfo.pos_z[i];
    j++;
    }
    for(k=0;k<live-1;k++){//live=2 or 3
      if((j-jstart)>(gridsize*2)) {printf("live=%d,j=%d,i=%d/%d,overflow\n",live,j,i,gridsize);exit(-1);}
      //else{
      x2[j-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
      y2[j-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
      z2[j-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
      j++;
      //}
    }
  }
  
  k = (j>banksize)*banksize + (j<=banksize)*j;
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+2*banksize+jstart,x2,sizeof(float)*(k-jstart), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+2*banksize+jstart,y2,sizeof(float)*(k-jstart), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+2*banksize+jstart,z2,sizeof(float)*(k-jstart), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaDeviceSynchronize());  
  free(sid1);
  free(x2);  free(y2);  free(z2);

  return j;
}
//==============================================================================
//=============== setbank() of __FTALLY2 samples to satisfy num_src ============
void setbank2(MemStruct DeviceMem, MemStruct HostMem, int banksize, unsigned jstart){
  jstart = (jstart>=banksize)*banksize + (jstart<banksize)*jstart;
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+banksize,DeviceMem.nInfo.pos_x+2*banksize,sizeof(float)*(jstart), cudaMemcpyDeviceToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+banksize,DeviceMem.nInfo.pos_y+2*banksize,sizeof(float)*(jstart), cudaMemcpyDeviceToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+banksize,DeviceMem.nInfo.pos_z+2*banksize,sizeof(float)*(jstart), cudaMemcpyDeviceToDevice));  
  
  if(jstart>=banksize) return;
  float* x2 = (float*)malloc(sizeof(float)*(banksize-jstart));
  float* y2 = (float*)malloc(sizeof(float)*(banksize-jstart));
  float* z2 = (float*)malloc(sizeof(float)*(banksize-jstart));

  int k; 
  for(k=jstart;k<banksize;k++){
    x2[k-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
    y2[k-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
    z2[k-jstart]=rand()*1.0/RAND_MAX*HostMem.wdspp[0];
  }

  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+banksize+jstart,x2,sizeof(float)*(banksize-jstart), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+banksize+jstart,y2,sizeof(float)*(banksize-jstart), cudaMemcpyHostToDevice));  
  gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+banksize+jstart,z2,sizeof(float)*(banksize-jstart), cudaMemcpyHostToDevice));  
  free(x2);  free(y2);  free(z2);
}

#endif
#endif


int count_pop(int *live, int gridsize){
  int sum = 0;
  for(int i=0;i<gridsize;i++)
    sum += (0!=live[i]);
  return sum;
}
#if defined(__1D)
void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin){
  int i=0;
  for(i=0;i<ubat;i++){//num_src is important as loop index, but useless in history<<<>>>
#if defined(__MTALLY)
    //printf("i=(%d/%d)",i,ubat);
#endif
    history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, num_src,i*gridx*blockx,banksize);
  }
  gpuErrchk(cudaDeviceSynchronize());  

#if defined(__MTALLY)
  //printf("\n");
#endif
}
#endif

#if defined(__3D)
//==============================================================================
//=========start_neutrons_active() launches kernels for active generations =====
//1. At least for __FTALLY2, source and delay bank can be performed segment-wise
void start_neutrons_active(unsigned ibat, unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_seg, unsigned banksize, unsigned tnum_bin, MemStruct HostMem){
  int i=0; 
  for(i=0;i<num_seg;i++){
    gpuErrchk(cudaMemset(DeviceMem.nInfo.imat, 0, gridx*blockx*sizeof(int)));
    history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, gridx*blockx*num_seg,i*gridx*blockx,banksize,1);
    gpuErrchk(cudaDeviceSynchronize()); 
    setbank_active_in(DeviceMem, HostMem, gridx*blockx,gridx*blockx*num_seg,i*gridx*blockx); 
    //printf("======After %d/%d segments, queue: ",i,num_seg);    HostMem.pull_list->ViewQueue();
  }
  setbank_active_out(DeviceMem, HostMem, gridx*blockx*num_seg);
  //printf("======After             , queue: ");     HostMem.pull_list->ViewQueue();
}


void     start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin){

  int i=0; 
  for(i=0;i<ubat;i++){//num_src is important as loop index, but useless in history<<<>>>
#if defined(__CTALLY2)
    gpuErrchk(cudaMemset(DeviceMem.cnt2_t, 0, tnum_bin*gridx*blockx*sizeof(int)));
#endif
    history<<<gridx, blockx/*, blockx*sizeof(unsigned)*/>>>(DeviceMem, num_src,i*gridx*blockx,banksize, 0);
    gpuErrchk(cudaDeviceSynchronize()); 
    //printf("%d/%d done\n",i,ubat);
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

#if defined(__CTALLY2)
  for(int i=0;i<num_bin;i++){
    reduce_sum_equal<<<gridx, blockx, blockx*sizeof(CMPTYPE)>>>(
                   DeviceMem.tally.cnt2+i*gridx*blockx, 
                   DeviceMem.block_spectrum+i*gridx);
  }
  for(int i=0;i<num_bin;i++){
    reduce_sum_equal<<<1, gridx, gridx*sizeof(CMPTYPE)>>>(
                   DeviceMem.block_spectrum+i*gridx, DeviceMem.batcnt2+i);
  }
  //printf("%s\n", cudaGetErrorString(cudaPeekAtLastError()));
  //printf("%s\n", cudaGetErrorString(cudaThreadSynchronize()));
  gpuErrchk(cudaMemcpy(HostMem.batcnt2,DeviceMem.batcnt2,sizeof(CMPTYPE)*num_bin, cudaMemcpyDeviceToHost));

#endif
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
