#ifndef __MATERIAL_H__
#define __MATERIAL_H__
#include "material_data.h"
#include "gpuerrchk.h"

class material{
public:
  unsigned numMat;
  unsigned *offsets;
  unsigned *isotopes;
  float* N_tot;
  float* densities;
public:
  material(struct matdata *, unsigned);
  ~material();
  void release_pointer();
};

#endif
