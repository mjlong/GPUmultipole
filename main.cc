#include "CPUComplex.h"
#include "hdf5IO.h"
#include "multipole_data.h"

extern void anyvalue(struct multipoledata, int*, double *);
int main(){
  char h5filename[] = "092238.h5";
  struct multipoledata U238data;
  int l;
  double d;
  h5read(U238data, h5filename);
  anyvalue(U238data,&l,&d);
  printf("numL=%d\n",l);
  printf("double=\%g\n",d);
  return 0;
}
