#ifndef __GPUERRCHK_H__
#define __GPUERRCHK_H__

#define gpuErrchk(ans){gpuAssert((ans),__FILE__, __LINE__);}
#include <cuda.h>
static void gpuAssert(cudaError_t code, char *file, int line, bool abort = true){
  if(cudaSuccess!=code){
    fprintf(stderr,"GPUassert:%s %s %d\n", cudaGetErrorString(code), file,line);
    if(abort) exit(code);
  }
}


#endif
