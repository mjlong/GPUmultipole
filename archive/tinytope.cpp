#include "tinytope.h"
isotope::isotope(char *filename){
  endfreadf2(filename);
}
isotope::~isotope(){
  free(pseudo_rho);
}
void isotope::endfreadf2(char* filename){
  int iL,iR,numr;
  double awri, factor,scattering_radius, channel_radius;
  FILE *file;
  char line[ENDFLEN];
  file = fopen(filename,"r");
  //Line 1 contains nothing relevant
  fgets(line, ENDFLEN, file);
  
  //Line 2 contains only ZAID and AWR
  awr = endfsci(line+11);
  
  //Line 3 contains nothing relevant
  fgets(line, ENDFLEN, file);
  
  //Line 4 contains E_low, E_high and NAPS
  fgets(line, ENDFLEN, file);
  flag_rcrs = endfint(line+55);
  
  //Line 5 contains target spin and number of L states
  fgets(line, ENDFLEN, file);
  number_l = endfint(line+44);
  //  allocate_l();
  pseudo_rho = (double*)malloc(number_l*sizeof(double));

  //Read in all the data
  for(iL=0;iL<number_l;iL++){
    fgets(line, ENDFLEN, file);
    awri=endfsci(line);
    scattering_radius  =endfsci(line+11);
    numr = endfint(line+55);
    //    initialize_l(iL);
    factor = awri/(awri+1.0);
    if((0.0!=scattering_radius)&&(1==flag_rcrs))
      channel_radius = scattering_radius;
    else
      channel_radius = (1.23*pow(awri,ONETRD) + 0.8)*0.1;
    pseudo_rho[iL]  = C1*factor*channel_radius;

    for(iR=0;iR<numr;iR++)
      fgets(line, ENDFLEN, file);
  }//end for iL

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

