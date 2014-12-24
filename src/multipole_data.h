#ifndef __MULTIPOLEDATA_H__
#define __MULTIPOLEDATA_H__

#define FILENAMELEN 20
#define MAXISOTOPES 10

#include <stdio.h>
#include <stdlib.h>
#include "CPUComplex.h"

#if defined(__CFLOAT)
#define CMPTYPE float
#else
#define CMPTYPE double
#endif

struct multipoledata{
  int fissionable;
  CPUComplex<CMPTYPE> *mpdata;
  unsigned      length;
  unsigned *l_value; // l index of the pole
  CMPTYPE   *pseudo_rho;  //inherit nomenclature from isotope.h

  int mode;                        // Spacing mode
  int windows;                     // Number of windows
  int fitorder;                    // Order of the fit. 1 linear, 2 quadratic, etc.
  int numL;                        // Number of l values
  CMPTYPE startE;                      // Start energy for the windows
  CMPTYPE endE;                        // End energy for the windows
  CMPTYPE spacing;                     // The actual spacing in the mode of choice.
  //Mode = 0 (linear)
  //spacing = inner
  //Mode = 1 (sqrt)
  //spacing = sqrt(multipole_w%endE - multipole_w%startE)/multipole_w%windows 
  //Mode = 2 (log)
  //spacing = log(multipole_w%endE - multipole_w%startE)/multipole_w%windows
  CMPTYPE sqrtAWR;
  int *w_start;// Contains the index of the pole at the start of the window
  int *w_end;  // Contains the index of the pole at the end of the window
  CMPTYPE *fit;
  //Contains the fitting function.  (reaction type, coeff index, window index)
  //=========================================================================
};

//after copy data from multipole_data to device, release memory
void freeMultipoleData(int numIsos, struct multipoledata* data);
unsigned count_isotopes(char* inputname);
void isotope_read(char* input, struct multipoledata* isotopes );

#if defined(__ALLCPU)
#define MP_EA 0
#define MP_RT 1
#define MP_RA 2
#define MP_RF 3
#define MP_FISS 1
#define FIT_T 0
#define FIT_A 1
#define FIT_F 2

#if defined(__MITW)
#include "Faddeeva.hh"
#define w_function Faddeeva::w
#endif
void host_xs_eval_fast(struct multipoledata iso, CMPTYPE E, CMPTYPE sqrtKT, 
			                 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);
void fill_factors(CMPTYPE sqrtE, int numL, CMPTYPE* pseudo_rho,   
                                        CPUComplex<double> *sigT_factor);
int pindex(int iP, int type);
int findex(int iW, int iC, int type, int orders, int types);
#endif

#endif
