#include "CPUComplex.h"
#include "hdf5IO.h"
#include "multipole_data.h"

extern void anyvalue(struct multipoledata, int*, double *, double *);
int main(){
  char h5filename[] = "092238.h5";
  struct multipoledata U238data;
  int l;
  double d1, d2;
  h5read(U238data, h5filename);
  anyvalue(U238data,&l,&d1, &d2);
  printf("numL=%d\n",l);
  printf("double=%10.6e and %10.6e\n",d1, d2);
  return 0;
}
