#include "isotope.h"
void isotope::endfreadf2(char* filename){
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
  int iL,iJ,iR;
  FILE *file;
  char line[ENDFLEN];
  resonance current_resonance;
  file = fopen(filename,"r");
  //Line 1 contains nothing relevant
  fgets(line, ENDFLEN, file);
  
  //Line 2 contains only ZAID and AWR
  fgets(line, ENDFLEN, file);
  ZAID = endfsci(line);
  //  atomic_weight_ratio = endfsci(line+11);
  
  //Line 3 contains nothing relevant
  fgets(line, ENDFLEN, file);
  
  //Line 4 contains E_low, E_high and NAPS
  fgets(line, ENDFLEN, file);
  E_low = endfsci(line);
  E_high = endfsci(line+11);
  flag_rcrs = endfint(line+55);
  
  //Line 5 contains target spin and number of L states
  fgets(line, ENDFLEN, file);
  target_spin = endfsci(line);
  //scattering_radius = endfsci(line+11);
  number_l = endfint(line+44);
  for(iL=0;iL<number_l;iL++){
    
  }
  printf("%d\n",number_l);
  fclose(file);
}

void isotope::initialize_l(){
  ljdegeneracy = (unsigned*)malloc(number_l*sizeof(unsigned)*2);
  scattering_radius = (double*)malloc(number_l*sizeof(double));
  atomic_weight_ratio = (double*)malloc(number_l*sizeof(double));
  channel_radius = (double*)malloc(number_l*sizeof(double));
  pseudo_k0 = (double*)malloc(number_l*sizeof(double));
  pseudo_k0r = (double*)malloc(number_l*sizeof(double));
  pseudo_k0r2 = (double*)malloc(number_l*sizeof(double));
  factor = (double*)malloc(number_l*sizeof(double));
  CONST = (double*)malloc(number_l*sizeof(double));
  pseudo_lambdabar2 = (double*)malloc(number_l*sizeof(double));
}
void isotope::check_degeneracy(){
  //Note: all agular momentum numbers are doubled to integer here
  unsigned *l_jdeg = (unsigned*)malloc(2*sizeof(unsigned));
  unsigned iL, iJ, offset=0;
  unsigned nspin = 1;
  unsigned tspin = (int)(2*target_spin);
  unsigned j_min;
  unsigned sum=0, deg;
  for(iL=0;iL<number_l;iL++){
    j_min = min(min(abs(2*iL+tspin-nspin),abs(2*iL-tspin+nspin)),abs(2*iL-tspin-nspin));
    l_jdeg[2*iL+1] = j_min; //odd indexes store j_min
    deg = iL + 0.5*nspin + target_spin - 0.5*j_min + 1; 
                     //even index stores deneneracy
    l_jdeg[2*iL] = deg;
    sum += deg;
  }
  degeneracy = (unsigned*)malloc(2*(number_l+sum)*sizeof(unsigned));
  //2*number_l  for each l, degeneracy and j_min are stored;
  //2*sum for each (l,j), number of channels and resonances are stored;
  degeneracy[iL]   = l_jdeg[2*iL];
  degeneracy[iL+1] = l_jdeg[2*iL+1];
  deg = l_jdeg[2*iL];
  for(iL=1;iL<number_l;iL++){
    offset += 2*deg+2;
    deg = l_jdeg[2*iL];
    degeneracy[offset]   = deg;
    degeneracy[offset+1] = l_jdeg[2*iL+1];
  }
  free(l_jdeg);
}
double endfsci(char *number){
  char sign = number[9];
  int i;
  int ten = number[10] - 48;
  int weithts[6]={100000,10000,1000,100,10,1};
  int value = 0;
  double decimal;
  for(i=0;i<6;i++)
    value += (number[i+3]-48)*weithts[i];
  decimal = number[1] - 48 + value*1.0E-6;
  if('+'==sign)
    for(i=1;i<=ten;i++)
      decimal *= 10.0;
  else
    for(i=1;i<=ten;i++)
      decimal *= 0.1;
  if('-' == number[0])
    decimal = 0.0 - decimal;
  return decimal;
}

unsigned endfint(char *number){
  int weights[10]={1E9,1E8,1E7,1E6,1E5,1E4,1000,100,10,1};
  int i, value=0;
  for (i=0;i<10;i++){
    if(' '!=number[i+1])
      value += (number[i+1]-48)*weights[i];
  }
  return value;
}

