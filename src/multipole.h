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
#define WINDOWS  4
#define DEVINTS  5

#define PLVAL 0
#define PPRHO 1
#define PWIND 2  //w_start and w_end share the same offset
#define PFITS 3  //fitT and fitA share the same offset
#define PFITF 4
#define PMPDATA 5
#define NUMOFFS 6

#define STARTE   0
#define SPACING  1
#define SQRTAWR  2
#define DEVREALS 4

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
#include "neutron.h"

#if defined(__MITW)
#include "Faddeeva.h"
#define w_function Faddeeva::w
#endif

#if defined(__QUICKW)
#include "QuickW.h"
#endif

#if defined(__FOURIERW)
#include "fourierw.h"
#endif
#include "global.h"
#include "gpuerrchk.h"
using namespace std;
//using namespace Faddeeva;

class multipole{
public:
  int *offsets;
  int *dev_numIso;
  int *dev_integers;
  CMPTYPE *dev_doubles;
  CComplex<CMPTYPE> *mpdata;
  unsigned *l_value; // l and j index of the pole
  //int w_function; //Which W function to use
  CMPTYPE  *pseudo_rho;  //inherit nomenclature from isotope.h

  int *w_start;// Contains the index of the pole at the start of the window
  int *w_end;  // Contains the index of the pole at the end of the window
  //CMPTYPE *fit;
  CMPTYPE *fitT;
  CMPTYPE *fitA;
  CMPTYPE *fitF;
  //Contains the fitting function.  (reaction type, coeff index, window index)
  //=========================================================================
#if defined(__QUICKWG) 
  CComplex<CMPTYPE>* mtable;
#endif
/*#if defined(__QUICKWC)
  CMPTYPE2 *mtable;
#endif*/
 public:
#if defined(__QUICKWG)
  multipole(struct multipoledata *data, int, CComplex<CMPTYPE> *wtable);
#else
  multipole(struct multipoledata *data, int);
#endif
  ~multipole();
  void release_pointer();
#if defined(__MITW) || defined(__QUICKW) || defined(__FOURIERW)
  __device__  void xs_eval_fast(int iM, CMPTYPE E, CMPTYPE sqrtKT, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);
#endif
#if defined(__SAMPLE)
  __device__  void xs_eval_fast(CMPTYPE E, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);
#endif
  /*__device__  void xs_eval_fast(CMPTYPE E, CMPTYPE sqrtKT, CMPTYPE rnd, 
					 CMPTYPE &sigT, CMPTYPE &sigA, CMPTYPE &sigF);*/
  __device__ void fill_factors(int prhoOffset, CMPTYPE sqrtE, int numL, CComplex<double> *sigT_factor);
                              //prhoOffset locates the pseudo_rho(iM,iL) in the long pseudo_rho array
  __host__ __device__  int findex(int, int, int, int, int);
  __host__ __device__  int pindex(int, int);

};
__device__ void broaden_n_polynomials(double En, double DOPP, double* factors, unsigned n);
#endif
