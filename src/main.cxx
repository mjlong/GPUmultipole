#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "process.h"

#include <time.h>
extern void createmptyh5(char *filename);
extern void writeh5_nxm_(char *filename, char *dsetname, double *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(char *filename, char *dsetname, int    *vec1, int *num_vec, int *length);
void printbless();
int main(int argc, char **argv){
  //printbless();
//============================================================ 
//====================calculation dimension===================
//============================================================
  unsigned print;
  unsigned gridx, blockx, gridsize,num_src;
  unsigned num_bin;
  unsigned num_bat;
  unsigned ubat,upto;
  double width, sigt, pf,pc;
  char name[50];
  int mode; //0=run only; 1=process only; 2=run & process
  if(argc>=8){//run or run+process
    gridx = atoi(argv[1]);
    blockx = atoi(argv[2]);
    gridsize = gridx*blockx;
    num_bat = atoi(argv[3]);    
    num_bin = atoi(argv[4]);
    width = atof(argv[5]);
    sigt  = atof(argv[6]);
    pf    = atof(argv[7]);
    pc    = atof(argv[8]);
    mode = 0;   //run only
    if(argc>=11){
      ubat = atoi(argv[9]);
      upto = atoi(argv[10]);
      print = atoi(argv[11]);
      mode = 2; //run+process
    }
  }
  else{
    mode = 1; //process only
    //name = argv[1];
    ubat = atoi(argv[2]);
    upto = atoi(argv[3]);
    print = atoi(argv[4]);
  }

  num_src=1;//num_src is not used but appears somewhere
//============================================================ 
//=============simulation memory allocation===================
//============================================================

  initialize_device();
  MemStruct HostMem, DeviceMem;

  initialize_memory(&DeviceMem, &HostMem, num_bin, gridx,blockx,num_bat,ubat);

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  copydata(DeviceMem,HostMem);
  printf("grid=[%3dx%3d],nhis=%-6d,ubat=%3d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridx,blockx,gridx*blockx,ubat,num_bat,num_bin,width);
  printf("mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));
//============================================================ 
//===============main simulation body=========================
//============================================================
  clock_t clock_start, clock_end;
  float time_elapsed = 0.f;
  unsigned active,banksize;
  active = 1;
  banksize = gridx*blockx;

  initialize_neutrons(gridx, blockx, DeviceMem,width); 
  clock_start = clock();

  for(int ibat=0;ibat<num_bat;ibat++){
  start_neutrons(gridx, blockx, DeviceMem, num_src,1,banksize);
  //active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);

  banksize = setbank(DeviceMem, gridsize);
  //printf("[%3d]%4d-->%4d: ", ibat,gridsize,banksize);
  save_results(ibat,gridx, blockx, num_src, num_bin, DeviceMem, HostMem);
  //resetcount(DeviceMem);
  resettally(DeviceMem.tally.cnt, num_bin*gridsize);
  }
  clock_end   = clock();
  time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
  printf("[time]  %d batches (*%d neutrons/batch) costs %f ms\n", num_bat,gridsize, time_elapsed);

  //============================================================================
  //==================== Write raw cnt to a hdf5 file ==========================
  //============================================================================
  char name1[10];  char name2[10];  char name3[10]; 
  sprintf(name1,"_%d",gridx*blockx);
  sprintf(name2,"_%d",ubat);
  sprintf(name3,"_%d",num_bat);
  int nbat=num_bat; int meshes=num_bin;
  strcpy(name,"Alltally"); strcat(name,name1); strcat(name,name2); strcat(name,name3); strcat(name,".h5");
  createmptyh5(name); //create empty file for future add dataset

  writeh5_nxm_(name,"batch_cnt", HostMem.batcnt, &nbat, &meshes);
  printf("[Save] Writing batch cnt to hdf5 .... done.\n");

  //============================================================================
  //=========================process the results ===============================
  //============================================================================

  clock_start = clock();
  //========================collison count to density ==========================
  cnt2flux(HostMem,gridsize,width/num_bin,num_bin,num_bat,ubat);
  //print_results(num_bin,num_bat,HostMem.acccnt);
  //print_results(num_bin,num_bat,HostMem.accmeans);
  printf("[Stat] Batch means and batch accmeans done:\n");
  if(0==print)
    print_results(num_bin,num_bat,HostMem.batchmeans);

  int intone=1;  int tempd;
  tempd = gridsize;  writeh5_nxm_(name,"num_history", &(tempd),  &intone, &intone);
  tempd = num_bat;   writeh5_nxm_(name,"num_batch",   &(tempd),  &intone, &intone);
  tempd = meshes;    writeh5_nxm_(name,"num_cells",   &(tempd),  &intone, &intone);
  writeh5_nxm_(name,"width",   &(width),  &intone, &intone);
  writeh5_nxm_(name,"sigma",   &(sigt),   &intone, &intone);
  writeh5_nxm_(name,"pf",      &(pf),     &intone, &intone);
  writeh5_nxm_(name,"pc",      &(pc),     &intone, &intone);

  //writeh5_nxm_(name,"batchmeans", HostMem.batchmeans, &nbat, &meshes);
  //printf("[Save] Writing means to hdf5 done:\n");

  //intone=nbat-ubat; writeh5_nxm_(name,"batchaccumu", HostMem.accmeans, &intone, &meshes);
  //printf("[Save] Writing acc means to hdf5 done:\n");

  //========================Average Square Error================================
  double *ASE = (double*)malloc(sizeof(double)*(num_bat-ubat));
  getASE(HostMem.accmeans, num_bin, num_bat,ubat, ref, ASE);
  printf("[Stat] ASE done:\n");
  if(0==print)
    print_results(num_bat-ubat,1,ASE);
  //=====================Auto-Correlation Coefficients==========================
  double *COR = (double*)malloc(sizeof(double)*upto*num_bin);
  getCOR(HostMem.batchmeans,num_bin,num_bat,ubat,upto,COR);
  printf("[Stat] Mesh correlations done:\n");
  if(0==print)
    print_results(upto,num_bin, COR);

  //==================== ACC fit ===============================================
  double *rho0s = (double*)malloc(sizeof(double)*num_bin);
  double *qs    = (double*)malloc(sizeof(double)*num_bin);
  fitall(COR,upto,num_bin,rho0s,qs);
  printf("[Stat] ACC fit done:\n");
  if(0==print){
    print_results(num_bin,1,rho0s);
    print_results(num_bin,1,qs);
  }
  //fitall1(COR,upto,num_bin,rho0s,qs);
  //printf("ACC fit done:\n");
  //print_results(num_bin,1,rho0s);
  //print_results(num_bin,1,qs);
  
  //=========================cell variance =====================================
  double *vars = (double*)malloc(sizeof(double)*num_bin);
  for(int im=0;im<num_bin;im++)
    vars[im] = variance(HostMem.batchmeans,num_bat,ubat,num_bin,im);
  printf("[Stat] Variance done:\n");
  if(0==print)
    print_results(num_bin,1,vars);


  //================= MASE (Mean average square error) =========================  
  double *EASE = (double*)malloc(sizeof(double)*(num_bat-ubat));
  getEASE(vars,num_bin,ubat,num_bat-ubat,rho0s,qs,EASE);
  printf("[Stat] EASE done:\n");
  if(0==print)
    print_results(num_bat-ubat,1,EASE);

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
  free(EASE);
  free(vars);
  free(rho0s);
  free(qs);
  free(COR);
  free(ASE);

  clock_end   = clock();
  time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
  printf("[time]  statistics costs %f ms\n", time_elapsed);


//============================================================ 
//=============simulation shut down===========================
//============================================================
  release_memory(DeviceMem, HostMem);
  printf("so far so good\n");

  return 0;
}







#define BLESS "[缪]"
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


