#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "process.h"

#include <time.h>
extern void createmptyh5(char *filename);
extern void createfixsrch5(char *filename);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, double *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, float  *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, int    *vec1, int *num_vec, int *length);
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
  int ubat,upto;
  int ibat=0;
  double width, sigt, pf,pc,v1;
  char name[50];
  int mode; //0=run only; 1=process only; 2=run & process
  int isSteady=0;

  unsigned gridr,blockr,gridsizr, ubatr;
  int num_batr,num_srcr;
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
  gridr = atoi(argv[10]);
  blockr = atoi(argv[11]);
  ubatr  = atoi(argv[12]);
  num_batr = atoi(argv[13]);


  if(0==gridsize){ 
    mode=0;     //read and run
  }
  else{
    mode=1;     //prepare and run
    createfixsrch5(argv[14]);
  }
  num_src=gridx*blockx*ubat;
  //gridr  = gridx/4;
  //blockr = blockx/4;
  gridsizr = blockr*gridr;
  //ubatr = ubat;
  num_srcr=gridr*blockr*ubatr;
  
  char name1[10];  char name2[10];  char name3[10]; 
  sprintf(name1,"_%d",gridsize);  sprintf(name2,"_%d",ubat);  sprintf(name3,"_%d",num_bat);
#if defined(__1D)
  strcpy(name,"R1dRawcnt"); 
#else
  strcpy(name,"R3dfixRawcnt"); 
#endif
  strcat(name,name1); strcat(name,name2); strcat(name,name3); 
  sprintf(name1,"_%d",gridsizr);  sprintf(name2,"_%d",ubatr);  sprintf(name3,"_%d",num_batr);
  strcat(name,name1); strcat(name,name2); strcat(name,name3); 
  strcat(name,".h5");
  createmptyh5(name); //create empty file for future add dataset
  
  int intone=1; 
  int inttwo=1;
  writeh5_nxm_(name,"/","num_batch_prep",  &(num_bat),  &intone, &intone);
  writeh5_nxm_(name,"/","num_batch",       &(num_batr),  &intone, &intone);
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
  initialize_memory(&DeviceMem, &HostMem, tnum_bin, gridx, blockx, ubat, gridr, blockr, ubatr);

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  HostMem.wdspp[5] = num_bin;
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  // note this only works for flat
  copydata(DeviceMem,HostMem);
  printf("nhis=%-6d,ubat=%3d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridsize,ubat,num_bat,num_bin,width);
  printf("mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));

//============================================================ 
//===============main simulation body=========================
//============================================================
  printf("[Info] Preparing fixed bank ... \n");
  unsigned active;
  int banksize;
  if(1==mode){//run fixed source preparation if fixed_source_file not specified
    active = 1;

    clock_start = clock();
    //==============================================================================
    //======================Steady State ===========================================
    //==============================================================================
    //======================Fixed source preparation ===============================
    banksize = gridx*blockx*ubat;
    initialize_neutrons(gridx, blockx, DeviceMem,width,banksize,ubat); 
    for(ibat=0;ibat<num_bat;ibat++){
      prep_neutrons(gridx, blockx, DeviceMem, ubat,1,banksize);
      //check(gridx,blockx,DeviceMem,ubat);
      //active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
      banksize = setbank(DeviceMem, HostMem, num_src);
      printf("[%3d]%4d-->%4d: \n", ibat,num_src,banksize);
    }
    clock_end   = clock();
    time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    printdone();
    printf("[time]  %d batches (*%d*%d neutrons/batch) costs %f ms\n", num_bat,gridsize,ubat, time_elapsed);
    writeh5_nxm_(argv[14],"/","gridsize",&gridsize,&intone,&intone);
    writeh5_nxm_(argv[14],"/","bunches", &ubat,    &intone,&intone);
    writeh5_nxm_(argv[14],"/","batches", &num_bat, &intone,&intone);
    writeh5_nxm_(argv[14],"/","x",       HostMem.nInfo.pos_x, &intone,&banksize);
    writeh5_nxm_(argv[14],"/","y",       HostMem.nInfo.pos_y, &intone,&banksize);
    writeh5_nxm_(argv[14],"/","z",       HostMem.nInfo.pos_z, &intone,&banksize);
  }
  else{//Reading from fixed_source_file
    

  }//end if (1==mode)     
    //====================== simulation with fixed source ======================
    clock_start = clock();
    printf("[Info] Running main simulation body ... \n");
    // plot initial distribution
    for(ibat=0;ibat<num_batr;ibat++){
      start_neutrons(gridr, blockr, DeviceMem, ubatr,num_src, banksize);
      //check(gridx,blockx,DeviceMem,ubat);
      printf("[%3d]%4d-->%4d: \n", ibat,banksize,num_srcr);
#if defined(__TALLY)
      save_results(ibat,gridr, blockr, tnum_bin, DeviceMem, HostMem);
      sprintf(name1,"%d",ibat);strcpy(name2,"batch_cnt");strcat(name2,name1);
      writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &tnum_bin);
      //resetcount(DeviceMem);
      resettally(DeviceMem.tally.cnt, tnum_bin*gridsizr);
#endif
      }

    clock_end   = clock();
    time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    printdone();
    printf("[time]  %d batches (*%d*%d neutrons/batch) costs %f ms\n", num_batr,gridsizr, ubatr,time_elapsed);

    //============================================================================
    //==================== Write raw cnt to a hdf5 file ==========================
    //============================================================================
    writeh5_nxm_(name, "/","num_history_prep",&(num_src) ,  &intone, &intone);
    writeh5_nxm_(name, "/","num_history",     &(num_srcr),  &intone, &intone);
    


  //============================================================================
  //=========================process the results ===============================
  //============================================================================

  /*

  FILE *fp=NULL;
  fp = fopen(name,"w");  
  fprintf(fp,"%.8e %.8e\n", gridx*blockx*1.0, num_bat*1.0-ubat);
  for(int i=0;i<num_bat-ubat;i++)
    fprintf(fp,"%.8e %.8e\n",ASE[i],EASE[i]);
  fclose(fp);
  

  strcpy(name,"boxtally");  strcat(name,name1);  strcat(name,name2);  strcat(name,name3);
  fp = NULL;
  fp = fopen(name,"w");  
  for(int i=0;i<num_bin;i++)
    fprintf(fp,"%.8e %.8e\n",HostMem.accmeans[(num_bat-ubat-1)*num_bin+i],ref);
  fclose(fp);

  //view ACC at cell print-1
  if(0<print){
    fp=NULL;
    strcpy(name,"acc");    sprintf(name1,"_%d",print-1);    strcat(name,name1);
    fp = fopen(name,"w");
    fprintf(fp,"%.8e\n",rho0s[print-1]);
    fprintf(fp,"%.8e\n",   qs[print-1]);
    for(int i=0;i<upto;i++)
      fprintf(fp,"%.8e\n",COR[(print-1)*upto+i]);
    fclose(fp);
  }

*/



//============================================================ 
//=============simulation shut down===========================
//============================================================
  release_memory(DeviceMem, HostMem);


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
  printf(" ..... done!\n");
}
