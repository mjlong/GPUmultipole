#ifndef __GPUERRCHK_H__
#define __GPUERRCHK_H__
#include <cuda.h>

#define gpuErrchk(ans){gpuAssert((ans),__FILE__, __LINE__);}
static void gpuAssert(cudaError_t code, char *file, int line, bool abort = true){
  if(cudaSuccess!=code){
    fprintf(stderr,"GPUassert:%s %s %d\n", cudaGetErrorString(code), file,line);
    if(abort) exit(code);
  }
}
#endif
/*
//issue of this error checkor:
//1. with asser(), error is reported with line in this file but not where it is called
#include <assert.h>
inline
cudaError_t gpuErrchk(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %sn", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}
*/

