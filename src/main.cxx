#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include <time.h>
extern void createfixsrch5(char *filename);
extern void createmptyh5(char *filename);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char *dsetname, double *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char *dsetname, float  *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char *dsetname, int    *vec1, int *num_vec, int *length);
extern void readh5_(char* filename, int* gridsize, int* nbat, 
	     int* meshes, double* width, 
	     double* sigt, double* pf, double* pc);
extern void readh5_(char* filename, int* cnt);
extern void readh5_(char* filename, float* x, float* y, float* z);

void printbless();
void printdone();
int main(int argc, char **argv){
  clock_t clock_start, clock_end;
  float time_elapsed = 0.f;
  //printbless();
//============================================================ 
//====================calculation dimension===================
//============================================================
  unsigned print;
  int gridx, blockx, gridsize,num_src;
  int num_bin, tnum_bin;
  int num_bat;
  int ubat;
  int ibat=0;
  double width, sigt, pf,pc,v1;
  char name[60];
  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
  num_bat = atoi(argv[3]);    
  num_bin = atoi(argv[4]);
  width = atof(argv[5]);
  sigt  = atof(argv[6]);
  pf    = atof(argv[7]);
  pc    = atof(argv[8]);
  ubat  = atoi(argv[9]);

  int num_ubat = num_bat%10000;
  num_bat = num_bat/10000;
  int delta_safe = atoi(argv[10])%1000;
  double Pf = pf/(pf+pc);
  printf("pf=%g,Pf=%g\n",pf,Pf);
  int delta_prep = (delta_safe); 
  int delta_dist = atoi(argv[11]);
  int num_seg = ubat;

  char name1[10];  char name2[10];  char name3[10]; char name4[10];
  sprintf(name1,"_%d",gridsize);
  sprintf(name2,"_%d",num_seg);
  sprintf(name3,"_%d",num_bat);
  sprintf(name4,"_%d",num_bin);
#if defined(__1D)&&defined(__FTALLY)
  strcpy(name,"R1dRawsrc"); 
#endif
#if defined(__1D)&&defined(__CTALLY)
  strcpy(name,"R1dRawcnt"); 
#endif

#if defined(__3D)&&(defined(__FTALLY)||defined(__FTALLY2))
  strcpy(name,"R3dUnRawsrc"); 
#endif
#if defined(__3D)&&defined(__CTALLY)
  strcpy(name,"R3d2Rawcnt_debug"); 
#endif
  strcat(name,name1); strcat(name,name2); strcat(name,name3); strcat(name,name4); 
  sprintf(name4,"_%dx%d_+%d_s%d",delta_dist,atoi(argv[10])%1000,num_ubat,atoi(argv[10])/1000);   strcat(name,name4); strcat(name,".h5");
  createmptyh5(name); //create empty file for future add dataset
  
  int intone=1; 
  int inttwo=1;
  writeh5_nxm_(name,"/","num_conv",   &(num_ubat),&intone, &intone);
  writeh5_nxm_(name,"/","num_prep", &(delta_prep),&intone, &intone);

  writeh5_nxm_(name,"/","num_batch",  &(num_bat),  &intone, &intone);
  writeh5_nxm_(name,"/","del_safe",&(delta_safe),  &intone, &intone);
  writeh5_nxm_(name,"/","num_cells",   &(num_bin),  &intone, &intone);
  writeh5_nxm_(name,"/","width",   &(width),  &intone, &intone);
  writeh5_nxm_(name,"/","sigma",   &(sigt),   &intone, &intone);
  writeh5_nxm_(name,"/","pf",      &(pf),     &intone, &intone);
  writeh5_nxm_(name,"/","pc",      &(pc),     &intone, &intone);
//============================================================ 
//=============simulation memory allocation===================
//============================================================
  //initialize_device();
  MemStruct HostMem, DeviceMem;
#if defined(__1D)
  tnum_bin = num_bin;
#endif
#if defined(__3D)
  tnum_bin = num_bin*num_bin*num_bin;
#endif
#if defined(__MTALLY)
  inttwo = tnum_bin*tnum_bin;
#endif
#if defined(__FTALLY)||defined(__FTALLY2)
  inttwo = tnum_bin;
#endif
  initialize_memory_data(&DeviceMem,&HostMem);
  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  HostMem.wdspp[5] = num_bin;
  HostMem.wdspp[6] = 3-2.0;  
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  // note this only works for flat
  copydata(DeviceMem,HostMem);
  printf("[]nhis=%-6d,num_seg=%d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridsize, num_seg,num_bat,num_bin,width);
  printf("[]mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));

//============================================================ 
//===============main simulation body=========================
//============================================================
  printf("[Info] Running main simulation body ... \n");
  int banksize;

  clock_start = clock();
    
  //==========================================================================
  //========== Converged source ==============================================
  //==========================================================================
  banksize = gridx*blockx*num_seg;
  num_src=gridx*blockx*num_seg;
  allocate_memory_converge(&DeviceMem, &HostMem, tnum_bin, gridx,blockx,num_seg);
  initialize_neutrons(gridx, blockx, DeviceMem,width,banksize,num_seg,atoi(argv[10])/1000); 

  for(ibat=0;ibat<num_ubat;ibat++){
    start_neutrons(gridx, blockx, DeviceMem, num_seg,num_src,banksize,tnum_bin);
    banksize=setbank_converge(DeviceMem, HostMem, num_src);
    //printf("%d[Converging source ...][%3d/%4d]%4d-->%4d: \n", -1,ibat,num_ubat,num_src,banksize);
  }
  //====================End of the phase to converge source: =====================
  //Note: the phase of convergence uses frame of size = gridsize*num_seg_XL, which is also used by the bank preparation phase
  //      and the last execution of setbank_convergence() has initialized the source for the next generation
  //      therefore, transition here is smooth
  char srcname[20];
  if( 0!=num_ubat ){
    sprintf(srcname,"srcpnts_%d",banksize);
    createfixsrch5(srcname);

    float* x2 = (float*)malloc(sizeof(float)*num_src*2);
    float* y2 = (float*)malloc(sizeof(float)*num_src*2);
    float* z2 = (float*)malloc(sizeof(float)*num_src*2);
    copysrcforwrite(HostMem, num_src, x2, y2, z2);
    writeh5_nxm_(srcname,"/","x",       x2, &intone,&banksize);
    writeh5_nxm_(srcname,"/","y",       y2, &intone,&banksize);
    writeh5_nxm_(srcname,"/","z",       z2, &intone,&banksize);
    writeh5_nxm_(srcname,"/","banksize",&banksize,           &intone,&intone);
    free(x2);  free(y2);  free(z2);

  }
  else{ //read
    int tempsize;readh5_(argv[12], &tempsize);
    printf("[Info] Reading fixed bank ...%d \n",tempsize);    
    float* x2 = (float*)malloc(sizeof(float)*tempsize);
    float* y2 = (float*)malloc(sizeof(float)*tempsize);
    float* z2 = (float*)malloc(sizeof(float)*tempsize);
    readh5_(argv[12], x2,y2,z2);
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+num_src,x2,sizeof(float)*tempsize, cudaMemcpyHostToDevice));  
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+num_src,y2,sizeof(float)*tempsize, cudaMemcpyHostToDevice));  
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+num_src,z2,sizeof(float)*tempsize, cudaMemcpyHostToDevice));  
    free(x2);  free(y2);  free(z2);
  }
  //==============================================================================
  //==========Start from the converged source; fill the delayed bank =============
  //==============================================================================
  unsigned delaysize = ((delta_prep+2)*gridx*blockx*num_seg);
  initialize_memory_bank(&HostMem, delaysize);
  //HostMem.pull_list = new CQueue<int>; 
  //HostMem.pull_list->InitQueue(num_src*2);
  HostMem.pull_list = new queue<int>; 

  for(ibat=0;ibat<delta_prep;ibat++){
    start_neutrons(gridx, blockx, DeviceMem, num_seg,num_src,banksize,tnum_bin);
    banksize=setbank_prepbank(DeviceMem, HostMem, num_src, ibat-delta_prep);
    printf("%d[Filling delay bank...][%3d/%4d]%4d-->%4d, cursor-->%d/%d: \n", -1,ibat,delta_prep,num_src,banksize, HostMem.bank.cursor_end[0],delaysize);
  }

  release_memory_converge(DeviceMem, HostMem);

  //==============================================================================
  //=========== Active generations ===============================================
  //==============================================================================

  //===================Reallocate memory in the frame of gridsize*num_seg ========
  banksize = gridx*blockx*num_seg;
  num_src=gridx*blockx*num_seg;
  allocate_memory_active(&DeviceMem, &HostMem, tnum_bin, gridx,blockx,num_seg);
  initialize_neutrons_active_not_src(gridx,blockx, DeviceMem,num_seg,atoi(argv[10])/1000+1);

  setbank_active_out(DeviceMem, HostMem, num_src);
  //printf("======After             , queue: ");  HostMem.pull_list->ViewQueue();

  clock_end   = clock();
  time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
  strcpy(name2,"timeprint0");
  writeh5_nxm_(name, "tally",name2, &(time_elapsed), &intone, &intone);

  for(ibat=0;ibat<num_bat;ibat++){
    clock_start = clock();
    start_neutrons_active(ibat, gridx, blockx, DeviceMem, num_seg,banksize,tnum_bin, HostMem);
    clock_end = clock();   time_elapsed += (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    sprintf(name1,"%d",ibat+1);strcpy(name2,"timeprint");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, &(time_elapsed), &intone, &intone);
    sprintf(name1,"%d",ibat);strcpy(name2,"BatcntA");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &inttwo);
    memset((HostMem).batcnt, 0, sizeof(CMPTYPE)*tnum_bin);
    printf("%d[Active tallying .....][%3d/%d]: \n", -1,ibat,num_bat);
  }

  //HostMem.pull_list->DeleQueue();
  delete HostMem.pull_list;
  release_memory_active(DeviceMem, HostMem);
  printdone();
  printf("[time]  %d batches (*%d neutrons/batch) costs %f ms\n", num_bat,gridsize, time_elapsed);

  //============================================================================
  //==================== Write raw cnt to a hdf5 file ==========================
  //============================================================================
  writeh5_nxm_(name, "/","num_history",&(num_src),  &intone, &intone);

//============================================================ 
//=============simulation shut down===========================
//============================================================

  release_memory_data(DeviceMem,HostMem);
  release_memory_bank(HostMem);

  return 0;
}







#define BLESS "[]"
void printbless(){
  printf(BLESS"                     _oo0oo_                  \n");  
  printf(BLESS"                    o8888888o                 \n");
  printf(BLESS"                    88\" . \"88                 \n");
  printf(BLESS"                    (| -_- |)                 \n"); 
  printf(BLESS"                    0\\  =  /0                 \n");
  printf(BLESS"                  ___/`---'\\___               \n");
  printf(BLESS"                .' \\\\|     |// '.             \n");
  printf(BLESS"               / \\\\|||  :  |||// \\            \n");
  printf(BLESS"              / _||||| -:- |||||- \\           \n");
  printf(BLESS"             |   | \\\\\\  -  /// |   |          \n");  
  printf(BLESS"             | \\_|  ''\\---/''  |_/ |          \n");
  printf(BLESS"             \\  .-\\__  '-'  ___/-. /          \n");
  printf(BLESS"           ___'. .'  /--.--\\  `. .'___        \n");
  printf(BLESS"        .\"\" \'<  `.___\\_<|>_/___.\' >\' \"\".      \n");
  printf(BLESS"       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |    \n");
  printf(BLESS"       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /    \n");
  printf(BLESS"   =====`-.____`.___ \\_____/___.-`___.-'===== \n");
  printf(BLESS"                     `=---='                  \n");
  printf(BLESS"\n\n\n");
  printf(BLESS"   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf(BLESS"        佛祖镇楼                  BUG辟易     \n");
}


void printdone(){
  printf("[..... done!]\n");
}
