#include "hdf5IO.h"
#include "multipole_data.h"

int main(){
  char h5filename[] = "092238.h5";
  struct multipoledata U238;
  h5read(U238, h5filename);
  printf("numL=%d\n", U238.length);
  return 0;
}
