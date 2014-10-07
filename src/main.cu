#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "global.h"
#include <stdio.h>
#include <string.h>

#define FILENAMELEN 20
#define MAXISOTOPES 10
/*#if defined (__QUICKWC)
#if defined(__CFLOAT)
__constant__ float2 table[LENGTH][LENGTH];
#else
__constant__ double2 table[LENGTH][LENGTH];
#endif
//__constant__ CMPTYPE table[LENGTH*LENGTH*2];
#endif
*/

#include <cuda.h>
#include "gpuerrchk.h"
extern void h5read(struct multipoledata & pole, char filename[]);
extern void anyvalue(struct multipoledata*,unsigned, unsigned, unsigned, unsigned, unsigned);
int init_data(char* input, char filenames[][FILENAMELEN]);

extern void tracemain(int num_particle, int, int, float*,float*, long long unsigned int);
void printbless();
int main(int argc, char **argv){
//test optix basic

  printbless();
  float geoPara[6] = {0.48f,0.5f,50.f,1.2f,100.f,100.f};
  //float geoPara[6] = {0.00048f,0.0005f,0.050f,0.0012f,0.100f,0.100f};
                      //r1,  r2,  h/2, p,   t,    H/2
  int argcs=5; 
  float *testmem; 
  CUdeviceptr my_ptr;
  unsigned width=atoi(argv[1+argcs]);
  gpuErrchk(cudaMalloc((void**)&testmem,sizeof(float)*width));
  gpuErrchk(cudaMalloc((void**)(&my_ptr), sizeof(float)*width));
  //cudaSetDeviceFlags(cudaDeviceMapHost|cudaDeviceLmemResizeToMax); 
  float testmemh[4]={7.2,3.2,4.1,5.9};
  gpuErrchk(cudaMemcpy(testmem, testmemh, sizeof(float)*width, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy((void*)my_ptr, testmemh, sizeof(float)*width, cudaMemcpyHostToDevice));
  tracemain(atoi(argv[1+argcs]), atoi(argv[2+argcs]), atoi(argv[3+argcs]), geoPara,testmem, my_ptr);
  gpuErrchk(cudaFree(testmem));

//end test optix

  int numIso;
  char filenames[MAXISOTOPES][FILENAMELEN];
  numIso = init_data(argv[5],filenames);
  struct multipoledata *isotopes;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  for(int i=0;i<numIso;i++)  
    h5read(isotopes[i],filenames[i]);

  anyvalue(isotopes,numIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));

  return 0;
}

int init_data(char* input, char line[][FILENAMELEN]){
  int numIso=-1;
  FILE *fp = NULL;
  fp = fopen(input,"r");

  if (fp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n",input);
    exit(1);
  }
  int i;
  while(fgets(line[++numIso], FILENAMELEN, fp) != NULL) {
    i=-1;
    printf("%s", line[numIso]);
    while('\n'!=line[numIso][++i]){}
    line[numIso][i]='\0';
    printf("\n"); 
  }
  fclose(fp);
  return numIso;
}









void printbless(){
  printf("                     _oo0oo_                  \n");  
  printf("                    o8888888o                 \n");
  printf("                    88\" . \"88                 \n");
  printf("                    (| -_- |)                 \n"); 
  printf("                    0\\  =  /0                 \n");
  printf("                  ___/`---'\\___               \n");
  printf("                .' \\\\|     |// '.             \n");
  printf("               / \\\\|||  :  |||// \\            \n");
  printf("              / _||||| -:- |||||- \\           \n");
  printf("             |   | \\\\\\  -  /// |   |          \n");  
  printf("             | \\_|  ''\\---/''  |_/ |          \n");
  printf("             \\  .-\\__  '-'  ___/-. /          \n");
  printf("           ___'. .'  /--.--\\  `. .'___        \n");
  printf("        ."" '<  `.___\\_<|>_/___.' >' "".      \n");
  printf("       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |    \n");
  printf("       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /    \n");
  printf("   =====`-.____`.___ \\_____/___.-`___.-'===== \n");
  printf("                     `=---='                  \n");
  printf("\n\n\n");
  printf("   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("        佛祖镇楼                  BUG辟易     \n");
}


