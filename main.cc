#include "CPUComplex.h"
#include "hdf5IO.h"
#include "multipole_data.h"

extern void anyvalue(struct multipoledata, int*, double *);
int main(){
  char h5filename[] = "092238.h5";
  struct multipoledata U238data;
  int l;
  double *ptr;
  h5read(U238data, h5filename);
  anyvalue(U238data,&l,ptr);
  printf("numL=%d\n",l);
  //  printf("pseudo_rhos\%g,%g\n",ptr[0],ptr[1]);
  return 0;
}
