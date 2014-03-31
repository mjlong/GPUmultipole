#ifndef __MULTIPOLE_H__
#define __MULTIPOLE_H__
#include "isotope.h"

#define MP_EA 1
#define MP_RT 2
#define MP_RA 3
#define MP_RF 4

class multipole{
 private:
  CComplex **mpdata; //All residues
  unsigned      npoles;
  unsigned *l_value, *j_value; // l and j index of the pole
  int      w_function; //Which W function to use
  double   atomic_weight_ratio, 
           startE, endE; //Start and End energy for windows
  double   *pseudo_rho;  //inherit nomenclature from isotope.h
  double   **gij;        //g statistic factor
  
 public:


};

#endif
