#include "h5_rdwt.h"
#include "multipole.h"
int main(){
  char h5filename[] = "092238.h5";
  double x,y,z;
  multipole U238;
  h5read(U238, h5filename);
  U238.xs_eval_fast(9.0,2.0,x,y,z);
  return 0;
}
