#include "CPUComplex.h"
#include "CComplex.h"
#include "multipole_data.h"
#include "multipole.h"

void anyvalue(struct multipoledata data, int *l, double *ptr){
  multipole U238(data);
  *l = U238.numL;
  //  ptr = U238.pseudo_rho;
}
