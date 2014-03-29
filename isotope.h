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
  /*========================================================================
    New Name            |  Old Name | Page #
    ZAID                |  ZA       | 49
    atomic_weight_ratio |  AWR      | 49
    E_low               |  EL       | 69
    E_high              |  EH       | 69
    number_channels     |  NC       | 52
    target_spin         |  SPI      | 72
    scattering_radius   |  AP       | 72
    number_l            |  NLS      | 72
    number_resonances   |  NRS      | 72
    resonance_energy    |  ER       | 72
    resonance_spin      |  AJ       | 72
    neutron_width       |  GN       | 73
    radiation_width     |  GG       | 73
    fission_width_1     |  GF       | 73
    fission_width_2
    scattering_radiuses |  APL      | 74
    flag_rcrs           | NASP      | 70
    ========================================================================*/
  unsigned    flag_rcrs, 
              number_l;
  double target_spin,         
         ZAID,                
         E_low,               
         E_high;
  //awrap stores l-dependent values: atomic_weight_ratio and scattering radius
  double *awrap; 
  resonance *resonances;
 public:
  void endfreadf2(char* filename);

};

unsigned endfint(char *number);
double endfsci(char *number);
#endif
