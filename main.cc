#include "CPUComplex.h"
#include "hdf5IO.h"
#include "multipole_data.h"

int main(){
  char h5filename[] = "092238.h5";
  struct multipoledata U238;
  h5read(U238, h5filename);
  U238.mpdata[0].output();
  U238.mpdata[1].output();
  return 0;
}
