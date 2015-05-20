#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "process.h"

#include <time.h>
extern void createmptyh5(char *filename);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, double *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, float  *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(char *filename, char *groupname, char *dsetname, int    *vec1, int *num_vec, int *length);
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
  int num_bin, tnum_bin;
  int num_bat;
  int ubat,upto;
  int ibat=0;
  double width, sigt, pf,pc,v1;
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
#if defined(__TRAN)
    v1 = atof(argv[9]);
    ubat = atoi(argv[10]);
#endif
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
#if defined(__1D)
  strcpy(name,"R1dRawcnt"); 
#else
#if defined(__TRAN)//3D transient
  strcpy(name,"R3tRawcnt"); 
#else              //3D steady state
  strcpy(name,"R3dRawcnt"); 
#endif
#endif
  strcat(name,name1); strcat(name,name3); strcat(name,".h5");
  createmptyh5(name); //create empty file for future add dataset
  
  int intone=1; 
  int inttwo=1;
  writeh5_nxm_(name,"/","num_batch",  &(num_bat),  &intone, &intone);
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
  initialize_memory(&DeviceMem, &HostMem, tnum_bin, gridx,blockx,num_bat,ubat);
#if defined(__PROCESS)
  if(1==mode)//process only, need to access the raw collision count
    readh5_(argv[1], HostMem.batcnt);
#endif    

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  HostMem.wdspp[5] = num_bin;
#if defined(__TRAN)
  CMPTYPE beta=0.01;
  CMPTYPE lambda = log(2.0)/10.0;
  CMPTYPE deltat = 1.0/(sigt*pf*2.5*v1);

  HostMem.wdspp[6] = v1;
  HostMem.wdspp[7] = beta;
  HostMem.wdspp[8] = lambda;
#endif
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  // note this only works for flat
  copydata(DeviceMem,HostMem);
  printf("nhis=%-6d,ubat=%3d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridsize,ubat,num_bat,num_bin,width);
  printf("mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));

//============================================================ 
//===============main simulation body=========================
//============================================================
  printf("[Info] Running main simulation body ... \n");
  unsigned active,banksize;
  if(1!=mode){//run simulation except 'process only' mode
    active = 1;

    clock_start = clock();
#if !defined(__TRAN)
    //==============================================================================
    //======================Steady State ===========================================
    //==============================================================================
    banksize = gridx*blockx;
    initialize_neutrons(gridx, blockx, DeviceMem,width,banksize); 
    // plot initial distribution
#if defined(__SCATTERPLOT)
    copyinitial(DeviceMem, HostMem, gridsize);
    sprintf(name1,"%d",ibat);
    strcpy(name2,"x");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
    strcpy(name2,"y");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
    strcpy(name2,"z");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
    strcpy(name2,"color");strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif
    for(ibat=0;ibat<num_bat;ibat++){
      start_neutrons(gridx, blockx, DeviceMem, num_src,1,banksize);
      //active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
      banksize = setbank(DeviceMem, HostMem, gridsize);
      printf("[%3d]%4d-->%4d: \n", ibat,gridsize,banksize);
#if defined(__TALLY)
      save_results(ibat,gridx, blockx, tnum_bin, DeviceMem, HostMem);
      sprintf(name1,"%d",ibat);strcpy(name2,"batch_cnt");strcat(name2,name1);
      writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &tnum_bin);
      //resetcount(DeviceMem);
      resettally(DeviceMem.tally.cnt, tnum_bin*gridsize);
#endif
#if defined(__SCATTERPLOT)
      sprintf(name1,"%d",ibat+1);
      strcpy(name2,"x");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
      strcpy(name2,"y");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
      strcpy(name2,"z");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
      strcpy(name2,"color");strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif
      }
#else 
    //==============================================================================
    //==================== Transient ===============================================
    //==============================================================================
    int allOld=0;
    banksize = gridsize;
    num_src = gridsize*ubat;
    printf("init....\n");
    initialize_neutrons(gridx, blockx, DeviceMem,width,banksize,num_src); 
    int csize = initialize_precursors(num_bat,banksize,num_src,lambda,beta*2.5*sigt*pf*v1,deltat,&HostMem);
    add_delayed_neutrons(DeviceMem, HostMem, 0, lambda, deltat, num_src,banksize);
//host_add_delayed(DeviceMem, HostMem, HostMem.initial_delayed[0], lambda, deltat, num_src,width);
//printf("%d addded\n",HostMem.initial_delayed[0]);
banksize += HostMem.initial_delayed[0];
    // plot initial distribution
#if defined(__SCATTERPLOT)
    copyinitial(DeviceMem, HostMem, gridsize);
    sprintf(name1,"%d",ibat);
    strcpy(name2,"live"); strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.live,   &intone, &gridsize);
    strcpy(name2,"x");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
    strcpy(name2,"y");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
    strcpy(name2,"z");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
    strcpy(name2,"color");strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif

    int *pops = (int*)malloc(sizeof(int)*num_bat);
    for(ibat=0;ibat<num_bat;ibat++){
      printf("ibat=%5d, banksize=%d\n",ibat,banksize);
      pops[ibat] = banksize;
      while(!allOld){
        transient_neutrons(gridx, blockx, DeviceMem, num_src,1,banksize,3-2.5*(1.0*(ibat<100)  + 1.002*(ibat>=100)));//don't worry, banksize is not used
#if defined(__TALLY)
	save_results(ibat,gridx,blockx, tnum_bin, DeviceMem, HostMem);
	sprintf(name1,"%d",ibat);strcpy(name2,"batch_cnt");strcat(name2,name1);
	writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &tnum_bin);
#endif
	banksize = flushbank(DeviceMem,HostMem,banksize,400.0,num_src,ibat,num_bat);
	//printf("flushed, banksize-->%d\n",banksize);
	allOld = (0==banksize);
      }
#if defined(__SCATTERPLOT)
      sprintf(name1,"%d",ibat+1);
      strcpy(name2,"live"); strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.live,   &intone, &gridsize);
      strcpy(name2,"x");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
      strcpy(name2,"y");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
      strcpy(name2,"z");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
      strcpy(name2,"color");strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif
      if(ibat<(num_bat-1)){
	add_delayed_neutrons(DeviceMem, HostMem, ibat+1, lambda, deltat, num_src,0);
	add_new_delayed(DeviceMem, HostMem, num_src, ibat, num_bat);
      }
      banksize = count_pop(HostMem.nInfo.live,num_src);

      allOld=0;
      if(0==banksize){
	printf("exiting: neutrons die out\n");
	break;
      }
#if defined(__TALLY)
      resettally(DeviceMem.tally.cnt, tnum_bin*num_src);
#endif
    }
#endif //end if transient

    clock_end   = clock();
    time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    printdone();
    printf("[time]  %d batches (*%d neutrons/batch) costs %f ms\n", num_bat,gridsize, time_elapsed);

    //============================================================================
    //==================== Write raw cnt to a hdf5 file ==========================
    //============================================================================
#if !defined(__TRAN)
    writeh5_nxm_(name, "/","num_history",&(gridsize),  &intone, &intone);
#else
    writeh5_nxm_(name, "/","neutron_pop",pops,  &intone, &num_bat);
    free(pops);
#endif

  }//end if (1!=mode) 

  //============================================================================
  //=========================process the results ===============================
  //============================================================================
#if defined(__PROCESS)
  if(0!=mode){//do process except 'run only' mode
    printf("[Info] Processing ... \n");
    strcpy(name,"RResult"); strcat(name,name1); strcat(name,name2); strcat(name,name3); strcat(name,".h5");
    createmptyh5(name); //create empty file for future add dataset
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
#endif
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
