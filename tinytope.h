#ifndef __TINYTOPE_H__
#define __TINYTOPE_H__ 

#include "stdio.h"
#include "stdlib.h"
#include <complex>
#include "global.h"

class isotope{
 private:
  unsigned    flag_rcrs, 
              number_l;
  double awr;
  /*========================================================================
    Calculated Values (Isotopic Only)
    
    factor - The ratio of the mass of the atom to the mass of the
             system of both the neutron and the atom.
    channel_radius - The calculated radius of the system.

    k0^2=2*mn *E/hbar^2              
    k ^2=2*mn'*E/hbar^2=k0^2*factor  (factor converts mn to reduced mass) 
    rho=k*a            (a = channel radius or scattering radius 
                             related with NAPS=1 or 0 in ENDF
                             (end102 page #321) )
    pseudo_k0    - k0 = pseuso_k0*sqrtE 
                 - pseudo_k0=sqrt(2mn/hb^2)[/pcm./sqrt(ev)]*factor
                 -          =0.0021968*factor[/pcm./sqrt(ev)]
    pseudo_rho0  - =pseudo_k0*channel_radius
    ========================================================================*/
  //l-dependent values: 
  double *pseudo_rho;

 public:
    isotope(){};
    isotope(char *filename);
    ~isotope();
    void endfreadf2(char* filename);
    friend class multipole;
};

unsigned endfint(char *number);
double endfsci(char *number);
#endif
