#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "global.h"

/*#if defined (__QUICKWC)
#if defined(__CFLOAT)
__constant__ float2 table[LENGTH][LENGTH];
#else
__constant__ double2 table[LENGTH][LENGTH];
#endif
//__constant__ CMPTYPE table[LENGTH*LENGTH*2];
#endif
*/
extern void h5read(struct multipoledata & pole, char filename[]);
extern void anyvalue(struct multipoledata*,unsigned, unsigned, unsigned, unsigned, unsigned);
int main(int argc, char *argv[]){
  unsigned numIso=2;
  char h5filename[] = "092238.h5";
  struct multipoledata *isotopes;//U238data;
  isotopes = (struct multipoledata*)malloc(sizeof(struct multipoledata)*numIso);
  for(int i=0;i<numIso;i++){
    h5read(isotopes[i], h5filename);
  }
  anyvalue(isotopes,numIso, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  return 0;
}
