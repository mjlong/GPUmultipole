#include "stdlib.h"
#include "CPUComplex.h"
#include "multipole_data.h"

extern void h5read(struct multipoledata & pole, char filename[]);
extern void anyvalue(struct multipoledata, int, int);
int main(int argc, char *argv[]){
  char h5filename[] = "092238.h5";
  struct multipoledata U238data;
  h5read(U238data, h5filename);
  anyvalue(U238data, atoi(argv[1]),atoi(argv[2]));
  //  printf("numL=%d\n",l);
  //  printf("double=%10.6e and %10.6e\n",d1, d2);
  return 0;
}
