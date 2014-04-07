#ifndef __MULTIPOLE_H__
#define __MULTIPOLE_H__

#include <complex>
using namespace std;
#define MP_EA 1
#define MP_RT 2
#define MP_RA 3
#define MP_RF 4

class multipole{
 private:
  int fissionable[1];
  complex<double> **mpdata; //All residues
  unsigned      length[1];
  unsigned *l_value, *j_value; // l and j index of the pole
  int      w_function; //Which W function to use
  double   atomic_weight_ratio; 
  double   *pseudo_rho;  //inherit nomenclature from isotope.h
  double   **gij;        //g statistic factor
  // Mode, set to 0 for linear, 1 for momentum, 2 for logarithmic.
  int mode[1];                        // Spacing mode
  int windows[1];                     // Number of windows
  int fitorder[1];                    // Order of the fit. 1 linear, 2 quadratic, etc.
  int numL[1];                        // Number of l values
  double startE[1];                      // Start energy for the windows
  double endE[1];                        // End energy for the windows
  double spacing[1];                     // The actual spacing in the mode of choice.
  //Mode = 0 (linear)
  //spacing = inner
  //Mode = 1 (sqrt)
  //spacing = sqrt(multipole_w%endE - multipole_w%startE)/multipole_w%windows 
  //Mode = 2 (log)
  //spacing = log(multipole_w%endE - multipole_w%startE)/multipole_w%windows
  int *w_start;// Contains the index of the pole at the start of the window
  int *w_end;  // Contains the index of the pole at the end of the window
  double *fit;
  //Contains the fitting function.  (reaction type, coeff index, window index)
    
  //=========================================================================
  // Storage Helpers
  int maxW;
  int ram;                         // RAM required for all 3 XS
  int ram_nofis;                   // RAM required for just 2


 public:
  multipole();
  multipole(char filename[]);
  void xs_eval_fast(double E, double sqrtKT, 
	       double *sigT, double *sigA, double *sigF);
  friend void h5read(multipole, char filename[]);
};

#endif
