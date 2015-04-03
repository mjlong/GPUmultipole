#include <stdio.h>
#include <string.h>

#include <neutron.h>
#include <manmemory.h>
#include "devicebridge.h"


#include <time.h>
void printbless();
int main(int argc, char **argv){
  //printbless();
//============================================================ 
//====================calculation dimension===================
//============================================================

  unsigned gridx, blockx, gridsize;
  unsigned num_src;

  gridx = atoi(argv[1]);
  blockx = atoi(argv[2]);
  gridsize = gridx*blockx;
  num_src = atoi(argv[3]);
  unsigned devstep;
  devstep = atoi(argv[4]);
//============================================================ 
//=============simulation memory allocation===================
//============================================================

  initialize_device();
  MemStruct HostMem, DeviceMem;
  unsigned num_bin = atoi(argv[6]);
  unsigned num_bat = atoi(argv[5]);
  initialize_memory(&DeviceMem, &HostMem, num_bin, gridx,blockx,num_bat);
  float width = atof(argv[7]);

  HostMem.wdspp[0] = width;
  HostMem.wdspp[1] = width/num_bin;
  HostMem.wdspp[2] = 1.0/atof(argv[8]); //sigmat
  HostMem.wdspp[3] = atof(argv[9]); //pf
  HostMem.wdspp[4] = atof(argv[10]);//pc
  copydata(DeviceMem,HostMem);
  printf("grid=[%3dx%3d],devstep=%3d,nhis=%-6d,nbat=%-6d,meshes=%-6d,box width=%.2f\n",gridx,blockx,devstep,num_src,num_bat,num_bin,width);
  printf("mfp=%.5f, pf=%.5f, pc=%.5f, ps=%.5f\n",HostMem.wdspp[2], HostMem.wdspp[3], HostMem.wdspp[4],1-(HostMem.wdspp[3]+HostMem.wdspp[4]));
//============================================================ 
//===============main simulation body=========================
//============================================================
  clock_t clock_start, clock_end;
  float time_elapsed = 0.f;
  unsigned active;
  active = 1;


  initialize_neutrons(gridx, blockx, DeviceMem,width); 
  clock_start = clock();
  start_neutrons(gridx, blockx, DeviceMem, num_src,1,devstep);
  active = count_neutrons(gridx, blockx, DeviceMem, HostMem,num_src);
 

  clock_end   = clock();
  time_elapsed = (float)(clock_end-clock_start)/CLOCKS_PER_SEC*1000.f;


  printf("[time], this batch costs %f ms\/%d neutrons\n", time_elapsed, HostMem.num_terminated_neutrons[0]);

  HostMem.num_terminated_neutrons[0]+=count_lives(gridx,blockx,DeviceMem,HostMem);

  print_results(0,gridx, blockx, num_src, num_bin, DeviceMem, HostMem, time_elapsed);

 

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


