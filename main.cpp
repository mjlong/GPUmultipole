#include "h5_rdwt.h"
#include "multipole.h"
int main(){
  char h5filename[] = "092238.h5";
  char endfname[]   = "endftest";
  double x,y,z;
  multipole U238;//(endfname);
  isotope U238i(endfname);
  U238.isotopeinfo(U238i);
  h5read(U238, h5filename);
  U238.xs_eval_fast(9.0,2.0,x,y,z);
  return 0;
}
