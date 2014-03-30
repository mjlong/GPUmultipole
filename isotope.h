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
   sqrtPF              - square root of the penetration factor
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
               sqrtPF,              
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
  //  double **numik;
  //TODO: numik seems only applies to Reich-Moore
  //========================================================================
  // For f-wave resonances only, contains the energy dependence    
  CComplex *QC; 
  CComplex QPF;
};




class isotope{
 private:
  bool NOSHIFT;
  unsigned    flag_rcrs, 
              number_l;
  double target_spin,         
         ZAID,                
         E_low,               
         E_high;
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
    pseudo_rho02 - =pseudo_rho0^2
     
    pseudo_lambdabar2 = pi^2/k^2*E
                      = pi^2*hb^2/(2mn)[b.ev]/factor^2
                      = 1301954.389/factor^2[b.ev]
    pseudo_twolambdabar2 = 2*pseudo_lambdabar2
    gij                  = (2J+1)/(2I+1)(2s+1) = 0.5*(2J+1)/(2I+1)
                         where s=neutron_spin = 0.5;
                         gij is statistic factor describing 
			 (2J+1) states out of 
			 (2l+1)(2I+1)(2s+1) states, 
			 (2l+1) disappears due to Legendre Polynomial
    ========================================================================*/
  //l-dependent values: 
  unsigned *l_jdeg;
  double *scattering_radius;
  double *atomic_weight_ratio;
  double *channel_radius;
  double *pseudo_k0;
  double *pseudo_rho0;
  double *pseudo_rho02;
  double *factor;
  double *pseudo_lambdabar2, *pseudo_twolambdabar2;
  // lj-dependent values:
  unsigned *number_channels, *number_resonances;
  double   *channel_spin, *gij;
  //ljr-dependent values
  resonance **resonances;
  /*========================================================================
    Calculated Values (Results)
    
    poles           - Positions of the poles in complex space
    residue_total   - Residues at the poles of the total XS
    residue_absorb  - Residues at the poles of the absorption XS
    residue_fission - Residues at the poles of the fission XS
    LJM             - TO-DO: Figure out exactly what this is.
    ========================================================================*/
    CComplex *poles,
             *residue_total,   
             *residue_absorb,  
             *residue_fission, 
             *LJM;
 public:
    isotope(){NOSHIFT=false;};
    void endfreadf2(char* filename);
    int check_degeneracy();
    void allocate_l();
    void initialize_l(int iL);
    void allocate_lj(int sum);
    void initialize_lj(int iL, int iJ, int iLJ);
    int index(int iL, int iJ);
    int index(int iL, double j);
    void assign_resonance(int iL, resonance**res_l);
    void set_resonance(int iL, int iJ, int iR);
    resonance & get_resonance(int iL, int iJ, int iR);
};

unsigned endfint(char *number);
double endfsci(char *number);
#endif
