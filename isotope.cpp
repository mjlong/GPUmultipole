#include "isotope.h"
void isotope::endfreadf2(char* filename){
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
  awrap = (double*)malloc(2*number_l*sizeof(double));

  printf("%d\n",number_l);
  fclose(file);
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
