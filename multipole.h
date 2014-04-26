#ifndef __MULTIPOLE_H__
#define __MULTIPOLE_H__

#define MP_EA 0
#define MP_RT 1
#define MP_RA 2
#define MP_RF 3
#define MP_FISS 1
#define FIT_T 0
#define FIT_A 1
#define FIT_F 2

#include <stdlib.h>
#include <stdio.h>
#include "CPUComplex.h"
#include "multipole_data.h"
#include "CComplex.h"
#include <iostream>
#include "Faddeeva.h"
#include "global.h"
using namespace std;
//using namespace Faddeeva;

class multipole{
 private:
  int fissionable;
  CComplex *mpdata;
  unsigned      length;
  unsigned *l_value, *j_value; // l and j index of the pole
  int      w_function; //Which W function to use
  double   *pseudo_rho;  //inherit nomenclature from isotope.h

  CComplex *sigT_factor;
  // Mode, set to 0 for linear, 1 for momentum, 2 for logarithmic.
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
  CComplex *Z_array, *W_array;
  //=========================================================================



 public:
  __host__ __device__  multipole(struct multipoledata data);
  __host__ __device__  multipole(char filename[]);
  //  void isotopeinfo(isotope );
  __host__ __device__  void xs_eval_fast(double E, double sqrtAWR, 
					 double &sigT, double &sigA, double &sigF);
  __host__ __device__  void xs_eval_fast(double E, 
					 double &sigT, double &sigA, double &sigF);
  __host__ __device__  void fill_factors(double sqrtE, double *twophi);
  __host__ __device__  int findex(int, int, int);
  __host__ __device__  int pindex(int, int);
  friend void h5read(multipole&, char filename[]);
};

#endif
