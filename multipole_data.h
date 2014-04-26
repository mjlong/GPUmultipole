#ifndef __MULTIPOLEDATA_H__
#define __MULTIPOLEDATA_H__

#include "CPUComplex.h"

struct multipoledata{
  int fissionable;
  CPUComplex *mpdata;
  unsigned      length;
  unsigned *l_value, *j_value; // l and j index of the pole
  double   *pseudo_rho;  //inherit nomenclature from isotope.h

  int mode;                        // Spacing mode
  int windows;                     // Number of windows
  int fitorder;                    // Order of the fit. 1 linear, 2 quadratic, etc.
  int numL;                        // Number of l values
  double startE;                      // Start energy for the windows
  double endE;                        // End energy for the windows
  double spacing;                     // The actual spacing in the mode of choice.
  //Mode = 0 (linear)
  //spacing = inner
  //Mode = 1 (sqrt)
  //spacing = sqrt(multipole_w%endE - multipole_w%startE)/multipole_w%windows 
  //Mode = 2 (log)
  //spacing = log(multipole_w%endE - multipole_w%startE)/multipole_w%windows
  double sqrtAWR;
  int *w_start;// Contains the index of the pole at the start of the window
  int *w_end;  // Contains the index of the pole at the end of the window
  double *fit;
  //Contains the fitting function.  (reaction type, coeff index, window index)
  //=========================================================================
};

#endif
