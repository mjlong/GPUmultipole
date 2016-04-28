#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"

#include "process.h"

#include <time.h>
extern void createfixsrch5(char *filename);
extern void createmptyh5(char *filename);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char
			 *dsetname, double *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char
			 *dsetname, float  *vec1, int *num_vec, int *length);
extern void writeh5_nxm_(const char *filename, const char *groupname, const char
			 *dsetname, int    *vec1, int *num_vec, int *length);
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
  int num_seg,upto;
  int ibat=0;
  double width, sigt, pf,pc,v1;
  char name[60];
  int mode; //0=run only; 1=process only; 2=run & process
  int isSteady=0;

  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
  num_bat = atoi(argv[3]);    
  num_bin = atoi(argv[4]);
  width = atof(argv[5]);
  sigt  = atof(argv[6]);
  pf    = atof(argv[7]);
  pc    = atof(argv[8]);
  num_seg  = atoi(argv[9]);

  int num_ubat = num_bat%10000;
  num_bat = num_bat/10000;

  num_src=gridx*blockx*num_seg;
  char name1[10];  char name2[10];  char name3[10]; char name4[10];
  sprintf(name1,"_%d",gridsize);
  sprintf(name2,"_%d",num_seg);
  sprintf(name3,"_%d",num_bat);
  sprintf(name4,"_%d",num_bin);
#if defined(__1D)&&defined(__MTALLY)
  strcpy(name,"R1d_UN_Tmacnt"); 
#endif
#if defined(__1D)&&defined(__FTALLY)
  strcpy(name,"R1dRawsrc"); 
#endif
#if defined(__1D)&&defined(__CTALLY)
  strcpy(name,"R1dRawcnt"); 
#endif
#if defined(__3D)&&defined(__MTALLY)
  strcpy(name,"R3d_UN_Tmacnt"); 
#endif
#if defined(__3D)&&(defined(__FTALLY)||defined(__FTALLY2))
  strcpy(name,"R3dRawsrc"); 
#endif
#if defined(__3D)&&defined(__CTALLY)
  strcpy(name,"R3d2Rawcnt_debug"); 
#endif
  strcat(name,name1); strcat(name,name2); strcat(name,name3); strcat(name,name4); 

#if defined(__MTALLY)
  sprintf(name4,"_s%d",atoi(argv[10]));
#else
  sprintf(name4,"_i%d_s%d",num_ubat,atoi(argv[10])); 
#endif
  
  strcat(name,name4); strcat(name,".h5");
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
#if defined(__MTALLY)
  inttwo = tnum_bin*tnum_bin;
#endif
#if defined(__FTALLY)||defined(__FTALLY2)
  inttwo = tnum_bin;
#endif
  initialize_memory(&DeviceMem, &HostMem, tnum_bin, gridx,blockx,num_bat,num_seg);

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/sigt;
  HostMem.wdspp[3] = pf;
  HostMem.wdspp[4] = pc;
  HostMem.wdspp[5] = num_bin;
  HostMem.wdspp[6] = 3-2.45;  
  double ref = 1.0/(HostMem.wdspp[3]+HostMem.wdspp[4])/width;
  // note this only works for flat
  copydata(DeviceMem,HostMem);
  printf("[]nhis=%-6d,num_seg=%3d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",
	 gridsize,num_seg,num_bat,num_bin,width);
  printf("[]mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));

//============================================================ 
//===============main simulation body=========================
//============================================================
  printf("[Info] Running main simulation body ... \n");
  int banksize, oldbanksize;
  clock_start = clock();
  //==============================================================================
  //======================Steady State ===========================================
  //==============================================================================
  banksize = gridx*blockx*num_seg;
#if defined(__MTALLY)||(__FTALLY_UN)
  banksize = banksize/2;
#endif
  initialize_neutrons(gridx, blockx, DeviceMem,width,banksize,num_seg,
		      atoi(argv[10])); 
  // plot initial distribution
#if defined(__SCATTERPLOT)
  copyinitial(DeviceMem, HostMem, gridsize);
  sprintf(name1,"%d",ibat);
  strcpy(name2,"x");    strcat(name2,name1); 
  writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
  strcpy(name2,"y");    strcat(name2,name1); 
  writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
  strcpy(name2,"z");    strcat(name2,name1); 
  writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
  strcpy(name2,"color");strcat(name2,name1); 
  writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif

#if !defined(__MTALLY)&&!defined(__FTALLY_UN)
  //==========================================================================
  //========== Converged source ==============================================
  //==========================================================================
  for(ibat=0;ibat<num_ubat;ibat++){
    start_neutrons(gridx, blockx, DeviceMem, num_seg,num_src,banksize,tnum_bin);
    banksize=setbank_converge(DeviceMem, HostMem, num_src);
    printf("[Converging source ...][%3d/%4d]%4d-->%4d: \n",
	   ibat,num_ubat,num_src,banksize);
  }

  char srcname[20];
  if( 0!=num_ubat ){
    sprintf(srcname,"srcpntM_%d",banksize);
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
    int tempsize;readh5_(argv[11], &tempsize);
    printf("[Info] Reading fixed bank ...%d \n",tempsize);    
    float* x2 = (float*)malloc(sizeof(float)*tempsize);
    float* y2 = (float*)malloc(sizeof(float)*tempsize);
    float* z2 = (float*)malloc(sizeof(float)*tempsize);
    readh5_(argv[11], x2,y2,z2);
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_x+num_src,x2,sizeof(float)
			 *tempsize, cudaMemcpyHostToDevice));  
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_y+num_src,y2,sizeof(float)
			 *tempsize, cudaMemcpyHostToDevice));  
    gpuErrchk(cudaMemcpy(DeviceMem.nInfo.pos_z+num_src,z2,sizeof(float)
			 *tempsize, cudaMemcpyHostToDevice));  
    free(x2);  free(y2);  free(z2);
  }

  //====================End of the phase to converge source ===================
#endif //not defined __MTALLY  and not defined __FTALLY
  clock_end   = clock();
  time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
  strcpy(name2,"timeprint0");
  writeh5_nxm_(name, "tally",name2, &(time_elapsed), &intone, &intone);
  for(ibat=0;ibat<num_bat;ibat++){
    oldbanksize = banksize;
    clock_start = clock();
#if (!defined(__FTALLY2))
    start_neutrons(gridx, blockx, DeviceMem, num_seg,num_src,banksize,tnum_bin);
#endif
    //check(gridx,blockx,DeviceMem,num_seg);
    //active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
#if defined(__TALLY)
#if defined(__MTALLY)||(__FTALLY)||(__FTALLY_UN)||(__FTALLY2)
#if !defined(__FTALLY2)
#if defined(__MTALLY)||(__FTALLY_UN)
    banksize = setbank(DeviceMem, HostMem, num_src,oldbanksize,tnum_bin);
#else
    banksize = setbank(DeviceMem, HostMem, num_src,tnum_bin);
#endif//end MTALLY
#else //else = defined FTALLY2
    banksize = start_neutrons(gridx, blockx, DeviceMem, num_seg,num_src,banksize,tnum_bin,HostMem);
#endif
    clock_end = clock();   time_elapsed += (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;
    sprintf(name1,"%d",ibat+1);strcpy(name2,"timeprint");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, &(time_elapsed), &intone, &intone);

    sprintf(name1,"%d",ibat); strcpy(name2,"Tmatrix");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &inttwo);
#if defined(__MTALLY)||(__FTALLY_UN)
    sprintf(name1,"%d",ibat);strcpy(name2,"sizeprint");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, &(banksize), &intone, &intone);
#endif

#if defined(__MTALLY)
    memset((HostMem).batcnt, 0, sizeof(CMPTYPE)*tnum_bin*tnum_bin);
#else
    memset((HostMem).batcnt, 0, sizeof(CMPTYPE)*tnum_bin);
#endif //end M or F
#else  //else C
    banksize = setbank(DeviceMem, HostMem, num_src);
    //initialize_neutrons_fix(gridx,blockx,DeviceMem,width,num_seg); 
    save_results(ibat,gridx, blockx, tnum_bin, DeviceMem, HostMem);
    sprintf(name1,"%d",ibat);strcpy(name2,"batch_cnt");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, HostMem.batcnt, &intone, &tnum_bin);
    resettally(DeviceMem.tally.cnt, tnum_bin*gridsize);
#if defined(__CTALLY2)
    sprintf(name1,"%d",ibat);strcpy(name2,"batch_cnts");strcat(name2,name1);
    writeh5_nxm_(name, "tally",name2, HostMem.batcnt2, &intone, &tnum_bin);
    resettally(DeviceMem.tally.cnt2, tnum_bin*gridsize);
#endif
#endif //end  TALLY types
    printf("%d[%3d]%4d-->%4d: \n", -1,ibat,oldbanksize,banksize);
#endif //end TALLY
#if defined(__SCATTERPLOT)
    sprintf(name1,"%d",ibat+1);
    strcpy(name2,"x");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_x,  &intone, &gridsize);
    strcpy(name2,"y");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_y,  &intone, &gridsize);
    strcpy(name2,"z");    strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.pos_z,  &intone, &gridsize);
    strcpy(name2,"color");strcat(name2,name1); writeh5_nxm_(name, "scatterplot",name2,HostMem.nInfo.energy, &intone, &gridsize);
#endif
  }

    printdone();
    printf("[time]  %d batches (*%d neutrons/batch) costs %f ms\n", num_bat,gridsize, time_elapsed);

    //============================================================================
    //==================== Write raw cnt to a hdf5 file ==========================
    //============================================================================
    writeh5_nxm_(name, "/","num_history",&(num_src),  &intone, &intone);


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
  printf("[..... done!]\n");
}
