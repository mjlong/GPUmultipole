#include "isotope.h"
void isotope::endfreadf2(char* filename){
  FILE *file;
  char line[ENDFLEN];
  resonance current_resonance;
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
  file = fopen(filename,"r");
  //Line 1 contains nothing relevant
  fgets(line, ENDFLEN, file);
  //Line 2 contains only ZAID and AWR
  fgets(line, ENDFLEN, file);
  
}

double endfsci(char number[11]){
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
  return decimal;
}
