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

#define MODE     0
#define FITORDER 1
#define NUML     2
#define FISSIONABLE  3
#define STARTE   0
#define SPACING  1
#define SQRTAWR  2

#if defined(__CFLOAT)
#define CMPTYPE float
#else
#define CMPTYPE double 
#endif
/*
  Mode, set to 0 for linear, 1 for momentum, 2 for logarithmic.
  Mode = 0 (linear)
  spacing = inner
  Mode = 1 (sqrt)
  spacing = sqrt(multipole_w%endE - multipole_w%startE)/multipole_w%windows 
  Mode = 2 (log)
  spacing = log(multipole_w%endE - multipole_w%startE)/multipole_w%windows
*/

#include <stdlib.h>
#include <stdio.h>
#include "CPUComplex.h"
#include "multipole_data.h"
#include "CComplex.h"
#include <iostream>

#if defined(__MITW)
#include "Faddeeva.h"
#endif

#if defined(__QUICKW)
#include "QuickW.h"
#endif

#include "global.h"
#include "gpuerrchk.h"
using namespace std;
//using namespace Faddeeva;

class multipole{
public:
  int *dev_integers;
  CMPTYPE *dev_doubles;
  CComplex<CMPTYPE> *mpdata;
  unsigned *l_value; // l and j index of the pole
  //int w_function; //Which W function to use
  CMPTYPE  *pseudo_rho;  //inherit nomenclature from isotope.h

  int *w_start;// Contains the index of the pole at the start of the window
  int *w_end;  // Contains the index of the pole at the end of the window
  CMPTYPE *fit;
  //Contains the fitting function.  (reaction type, coeff index, window index)
  //=========================================================================
#if defined(__QUICKW)
  CComplex<CMPTYPE>* table;
#endif

 public:
#if defined(__QUICKW)
  multipole(struct multipoledata data, CComplex<CMPTYPE> *wtable);
#else
  multipole(struct multipoledata data);
#endif
  ~multipole();
  void release_pointer();
#if defined(__MITW) || defined(__QUICKW)
  __device__  void xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);
#endif
  __device__  void xs_eval_fast(CMPTYPE E, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);
  /*__device__  void xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, CMPTYPE rnd, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);*/
  __device__ void fill_factors(CMPTYPE sqrtE, int numL, CComplex<CMPTYPE> *sigT_factor);
  __host__ __device__  int findex(int, int, int, int, int);
  __host__ __device__  int pindex(int, int);

};

#endif
