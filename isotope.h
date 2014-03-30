#ifndef __ISOTOPE_H__
#define __ISOTOPE_H__ 

#include "stdio.h"
#include "stdlib.h"
#include "CPUComplex.h" //TODO:update to CComplex.h
#include "global.h"
struct resonance{
  /*======================================================================
   From ENDF
   
   E               - Energy of the center of the resonance
   sqrtE           - Square root of the energy 
   radiation_width - Radiation width
   B     - Array containing, in order:
             Neutron Width
             Fission Width #1
             Fission Width #2
  ======================================================================*/
  
  double E, radiation_width, neutron_width, fission_width_1, fission_width_2;
  /*========================================================================
   Calculated Values (Multipole)
   
   rho   - pseudo_k0*sqrtE*channel_radius
   rho2  - rho squared
   rho4  - rho to the fourth
   rho6  - rho to the sixth (it gets used frequently, may as well
           calculate it once and be done with it)
  
   penetration_factor - The penetration factor
   sqtPF              - square root of the penetration factor
   shift_factor       - The shift factor.
   ER                 - E - HALFI*gamma
   neutron_width_0    - neutron_width/sqrtE, or the neutron width
                        that is at the energy E.
   fission_width      - the sum of the fission widths
   total_width        - the sum of all widths at energy E
   absorption_width   - same, but only for radiation/fission
   EA                 - E - HALFI*total_width
   COE1               - -HALFI*neutron_width/sqrtE
   COE2               - E - HALFI*absorption_width
   COEF2              - TO-DO
  ======================================================================*/
  double       rho,                
               rho2,               
               rho4,               
               rho6,               
               sqrtE,              
               penetration_factor, 
               sqtPF,              
               shift_factor,       
               fission_width,      
               total_width,        
               absorption_width,   
               neutron_width_0,    
               B[3] ;
  CComplex        ER,                  
                  EA,              
                  COE1,            
                  COE2,            
                  COEF2;
  //========================================================================
  //NUMIK( Numerator of the last term of eq D14 (I-K) from ENDF-102 )
  double **numik;
  //========================================================================
  // For f-wave resonances only, contains the energy dependence    
  CComplex *QC; 
  CComplex QPF;
};
class isotope{
 private:
  unsigned    flag_rcrs, 
              number_l;
  double target_spin,         
         ZAID,                
         E_low,               
         E_high;
  //l-dependent values: 
  unsigned *l_jdeg;
  double *scattering_radius;
  double *atomic_weight_ratio;
  double *channel_radius;
  double *pseudo_k0;
  double *pseudo_k0r;
  double *pseudo_k0r2;
  double *factor;
  double *CONST;
  double *pseudo_lambdabar2;
  // lj-dependent values:
  //degeneracy stores jdegeneracy(l), jmin(l), number_resonances(l,j), number_channels(l,j)
  unsigned *number_channels, *number_resonances;
  double   *channel_spin, *gij;
  
  
  //ljr-dependent values
  resonance *resonances;
 public:
  void endfreadf2(char* filename);
  int check_degeneracy();
  void initialize_l();
  void initialize_lj(int sum);

};

unsigned endfint(char *number);
double endfsci(char *number);
#endif
