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
extern void readh5_(char* filename, int* gridsize, int* nbat, 
	     int* meshes, double* width, 
	     double* sigt, double* pf, double* pc);
extern void readh5_(char* filename, int* cnt);

void printbless();
void printdone();
int main(int argc, char **argv){
  clock_t clock_start, clock_end;
  float time_elapsed = 0.f;
  printbless();
//============================================================ 
//====================calculation dimension===================
//============================================================
  unsigned print;
  int gridx, blockx, gridsize,num_src;
  int num_bin;
  int num_bat;
  int ubat,upto;
  double width, sigt, pf,pc;
  char name[50];
  int mode; //0=run only; 1=process only; 2=run & process
  int isSteady=0;
  if(argc>=8+1){//run or run+process
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
    if(argc>=11+1){
      ubat = atoi(argv[9]);
      upto = atoi(argv[10]);
      print = atoi(argv[11]);
      mode = 2; //run+process
    }
  }
  else{
    mode = 1; //process only
    ubat = atoi(argv[2]);
    upto = atoi(argv[3]);
    print = atoi(argv[4]);
    readh5_(argv[1],&gridsize,&num_bat,&num_bin,&width,&sigt,&pf,&pc);
  }

  num_src=1;//num_src is not used but appears somewhere
  char name1[10];  char name2[10];  char name3[10]; 
  sprintf(name1,"_%d",gridsize);
  sprintf(name2,"_%d",ubat);
  sprintf(name3,"_%d",num_bat);


//============================================================ 
//=============simulation memory allocation===================
//============================================================
  initialize_device();
  MemStruct HostMem, DeviceMem;
#if defined(__1D)
  initialize_memory(&DeviceMem, &HostMem, num_bin, gridx,blockx,num_bat,ubat);
#endif
#if defined(__3D)
  int tnum_bin = num_bin*num_bin*num_bin;
  initialize_memory(&DeviceMem, &HostMem, tnum_bin, gridx,blockx,num_bat,ubat);
#endif
  if(1==mode)//process only, need to access the raw collision count
    readh5_(argv[1], HostMem.batcnt);
    

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  // note this only works for flat
  copydata(DeviceMem,HostMem);
  printf("nhis=%-6d,ubat=%3d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridsize,ubat,num_bat,num_bin,width);
  printf("mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));
  int intone=1; 
  int inttwo=1;

//============================================================ 
//===============main simulation body=========================
//============================================================
  printf("[Info] Running main simulation body ... \n");
  unsigned active,banksize;
  if(1!=mode){//run simulation except 'process only' mode
    active = 1;

    clock_start = clock();
#if !defined(__TRAN)
    banksize = gridx*blockx;
    initialize_neutrons(gridx, blockx, DeviceMem,width,banksize); 
    for(int ibat=0;ibat<num_bat;ibat++){
      start_neutrons(gridx, blockx, DeviceMem, num_src,1,banksize);
      //active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);

      banksize = setbank(DeviceMem, gridsize);
      //printf("[%3d]%4d-->%4d: ", ibat,gridsize,banksize);
      save_results(ibat,gridx, blockx, num_bin, DeviceMem, HostMem);
      //resetcount(DeviceMem);
      resettally(DeviceMem.tally.cnt, num_bin*gridsize);
    }
#else
    int allOld=0;
    banksize = 40;
    initialize_neutrons(gridx, blockx, DeviceMem,width,banksize); 
    int *pops = (int*)malloc(sizeof(int)*num_bat);
    for(int ibat=0;ibat<num_bat;ibat++){
      pops[ibat] = banksize;
      while(!allOld){
	transient_neutrons(gridx, blockx, DeviceMem, num_src,1,banksize);
	save_results(ibat,gridx,blockx, tnum_bin, DeviceMem, HostMem);
	banksize = flushbank(DeviceMem,HostMem,banksize,400.0,gridsize);
	allOld = (0==banksize);
      }
      banksize = count_pop(HostMem.nInfo.live,gridsize);
      allOld=0;
      if(0==banksize){
	printf("exiting: neutrons die out\n");
	break;
      }
    }
#endif

    clock_end   = clock();
    time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    printdone();
    printf("[time]  %d batches (*%d neutrons/batch) costs %f ms\n", num_bat,gridsize, time_elapsed);

    //============================================================================
    //==================== Write raw cnt to a hdf5 file ==========================
    //============================================================================

    printf("[Save] Writing batch cnt to hdf5 .... ");
    strcpy(name,"RRawcnt"); strcat(name,name1); strcat(name,name3); strcat(name,".h5");
    createmptyh5(name); //create empty file for future add dataset
    writeh5_nxm_(name,"batch_cnt", HostMem.batcnt, &num_bat, &tnum_bin);
    printdone();
#if !defined(__TRAN)
    writeh5_nxm_(name,"num_history", &(gridsize),  &intone, &intone);
#else
    writeh5_nxm_(name,"neutron_pop", pops,  &intone, &num_bat);
    free(pops);
#endif
    writeh5_nxm_(name,"num_batch",   &(num_bat),  &intone, &intone);
    writeh5_nxm_(name,"num_cells",   &(num_bin),  &intone, &intone);
    writeh5_nxm_(name,"width",   &(width),  &intone, &intone);
    writeh5_nxm_(name,"sigma",   &(sigt),   &intone, &intone);
    writeh5_nxm_(name,"pf",      &(pf),     &intone, &intone);
    writeh5_nxm_(name,"pc",      &(pc),     &intone, &intone);

  }//end if (1!=mode) 

  //============================================================================
  //=========================process the results ===============================
  //============================================================================
  if(0!=mode){//do process except 'run only' mode
    printf("[Info] Processing ... \n");
    strcpy(name,"RResult"); strcat(name,name1); strcat(name,name2); strcat(name,name3); strcat(name,".h5");
    createmptyh5(name); //create empty file for future add dataset
    writeh5_nxm_(name,"num_history", &(gridsize),  &intone, &intone);
    writeh5_nxm_(name,"num_batch",   &(num_bat),  &intone, &intone);
    writeh5_nxm_(name,"num_cells",   &(num_bin),  &intone, &intone);
    writeh5_nxm_(name,"width",   &(width),  &intone, &intone);
    writeh5_nxm_(name,"sigma",   &(sigt),   &intone, &intone);
    writeh5_nxm_(name,"pf",      &(pf),     &intone, &intone);
    writeh5_nxm_(name,"pc",      &(pc),     &intone, &intone);
    writeh5_nxm_(name,"num_ubat",&(ubat),   &intone, &intone);
    writeh5_nxm_(name,"num_acc", &(upto),   &intone, &intone);

    clock_start = clock();
    //========================collison count to density ==========================
    printf("[Stat] Batch means and batch accmeans .... ");
    cnt2flux(HostMem,gridsize,width/num_bin,num_bin,num_bat,ubat);
    printdone();
    if(0==print)
      print_results(num_bin,num_bat,HostMem.batchmeans);
    //----------------------------------------------------------------------------
    printf("[Save] Writing means to hdf5... ");
    writeh5_nxm_(name,"batchmeans", HostMem.batchmeans, &num_bat, &num_bin);
    printf("... writing acc means to hdf5... ");
    intone=num_bat-ubat; writeh5_nxm_(name,"batchaccumu", HostMem.accmeans, &intone, &num_bin);
    printdone();
    //========================Average Square Error================================
    printf("[Stat] Average Square Error ... ");
    double *ASE = (double*)malloc(sizeof(double)*(num_bat-ubat));
    getASE(HostMem.accmeans, num_bin, num_bat,ubat, ref, ASE);
    printdone();
    if(0==print)
      print_results(num_bat-ubat,1,ASE);
    //----------------------------------------------------------------------------
    printf("[Save] Writing ASE to hdf5... ");
    inttwo=num_bat-ubat; intone=1; writeh5_nxm_(name,"ASE", ASE, &intone, &inttwo);
    printdone();
    //=====================Auto-Correlation Coefficients==========================
    printf("[Stat] Auto-correlation coefficients ... ");
    double *COR = (double*)malloc(sizeof(double)*upto*num_bin);
    getCOR(HostMem.batchmeans,num_bin,num_bat,ubat,upto,COR);
    printdone();
    if(0==print)
      print_results(upto,num_bin, COR);
    //----------------------------------------------------------------------------
    printf("[Save] Writing ACC to hdf5... ");
    writeh5_nxm_(name,"ACC", COR, &num_bin, &upto);
    printdone();

    //==================== ACC fit ===============================================
    printf("[Stat] ACC fit...");
    double *rho0s = (double*)malloc(sizeof(double)*num_bin);
    double *qs    = (double*)malloc(sizeof(double)*num_bin);
    fitall(COR,upto,num_bin,rho0s,qs);
    printdone();
    //fitall1(COR,upto,num_bin,rho0s,qs);
    //printf("ACC fit done:\n");
    //print_results(num_bin,1,rho0s);
    //print_results(num_bin,1,qs);
    if(0==print){
      print_results(num_bin,1,rho0s);
      print_results(num_bin,1,qs);
    }
    //----------------------------------------------------------------------------
    printf("[Save] Writing ACC fit result to hdf5... ");
    intone=1;
    writeh5_nxm_(name,"rho0s", rho0s, &intone, &num_bin);
    writeh5_nxm_(name,"qs",    rho0s, &intone, &num_bin);
    printdone();
  
    //=========================cell variance =====================================
    printf("[Stat] Variance ....");
    double *vars = (double*)malloc(sizeof(double)*num_bin);
    for(int im=0;im<num_bin;im++)
      vars[im] = variance(HostMem.batchmeans,num_bat,ubat,num_bin,im);
    printdone();
    if(0==print)
      print_results(num_bin,1,vars);
    //----------------------------------------------------------------------------
    printf("[Save] Writing mesh variances to hdf5... ");
    intone=1; writeh5_nxm_(name,"var", vars, &intone, &num_bin);
    printdone();

    //================= MASE (Mean average square error) =========================  
    printf("[Stat] Expected Average Square Error ....");
    double *EASE = (double*)malloc(sizeof(double)*(num_bat-ubat));
    getEASE(vars,num_bin,ubat,num_bat-ubat,rho0s,qs,EASE);
    printdone();
    if(0==print)
      print_results(num_bat-ubat,1,EASE);
    printf("[Save] Writing EASE to hdf5 ... ");
    intone=1; inttwo=num_bat-ubat; writeh5_nxm_(name,"EASE", EASE, &intone, &inttwo);
    printdone();

    free(EASE);
    free(vars);
    free(rho0s);
    free(qs);
    free(COR);
    free(ASE);

    clock_end   = clock();
    time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    printf("[time]  statistics costs %f ms\n", time_elapsed);

  }//end if(0!=mode) //end process
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


void printdone(){
  printf(" ..... done!\n");
}
