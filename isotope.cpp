#include "isotope.h"
isotope::isotope(char *filename){
  endfreadf2(filename);
}
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
  int iL,iJ,iLJ,iR,iJR,jdeg,sum,numr;
  double current_j, jmin;
  FILE *file;
  char line[ENDFLEN];
  resonance current_resonance;
  resonance **res_l;
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
  allocate_l();
  sum = check_degeneracy();
  allocate_lj(sum);

  //Read in all the data
  for(iL=0;iL<number_l;iL++){
    jdeg = l_jdeg[2*iL];
    fgets(line, ENDFLEN, file);
    atomic_weight_ratio[iL]=endfsci(line);
    scattering_radius[iL]  =endfsci(line+11);
    numr = endfint(line+55);
    initialize_l(iL);
    res_l=(resonance**)malloc(jdeg*sizeof(resonance*));
    for(iJ=0;iJ<l_jdeg[2*iL];iJ++){
      iLJ = index(iL,iJ);
      res_l[iJ]=(resonance*)malloc(numr*sizeof(resonance));
      number_channels[iLJ] = 1;
      number_resonances[iLJ] = 0;
    }//end for iJ
    //Read in resonances of iL
    jmin = l_jdeg[2*iL+1]*0.5;
    for(iR=0;iR<numr;iR++){
      fgets(line, ENDFLEN, file);
      current_resonance.E = endfsci(line);
      current_j = endfsci(line+11);
      iJ  = current_j - jmin ;
      iLJ = index(iL, iJ);
      initialize_lj(iL,iJ,iLJ);
      //iLJ depends on iL,iJ, but why not using as given above
      current_resonance.neutron_width = endfsci(line+22);
      current_resonance.radiation_width = endfsci(line+33);
      current_resonance.fission_width_1 = endfsci(line+44);
      current_resonance.fission_width_2 = endfsci(line+55);
      number_channels[iLJ] += ((0!=current_resonance.fission_width_1) +
			       (0!=current_resonance.fission_width_2));
      *(res_l[iJ]+number_resonances[iLJ]++) = current_resonance;
    }//end reading numr resonances of iL
    assign_resonance(iL, res_l);//res_l elements are freed
    free(res_l);                //res_l itself is freed
  }//end for iL
 

  fclose(file);
}
void isotope::assign_resonance(int iL, resonance** res_l){
  int jdeg = l_jdeg[2*iL];
  int ind0 = index(iL,(int)0);
  int ind1 = index(iL,jdeg-1);
  int iJ=0,numr;
  for(int iLJ=ind0;iLJ<=ind1;iLJ++){
    //ind0 adds to ind1 while (iJ=)0 adds to jdeg-1
    numr = number_resonances[iLJ];
    resonances[iLJ]=(resonance*)malloc(numr*sizeof(resonance));
    for(int iR=0;iR<numr;iR++){
      *(resonances[iLJ]+iR) = *(res_l[iJ]+iR);
      set_resonance(iL, iJ, iR);
    }
    free(res_l[iJ++]);
  }
}

void isotope::set_resonance(int iL, int iJ, int iR){
  resonance res;
  res = get_resonance(iL,iJ,iR);
  res.sqrtE = sqrt(abs(res.E));
  res.rho = pseudo_k[iL]*res.sqrtE*channel_radius[iL];
  res.rho2= res.rho * res.rho;
  res.rho4= res.rho2* res.rho2;
  res.rho6= res.rho2* res.rho4;
  //Set channel widths
  res.fission_width = abs(res.fission_width_1) + abs(res.fission_width_2);
  res.total_width = res.fission_width + res.neutron_width + res.radiation_width;
  res.absorption_width = res.fission_width + res.radiation_width;
  res.neutron_width_0 = res.neutron_width / res.sqrtE;
  res.B[0] = sqrt(abs(res.neutron_width/res.sqrtE))*((int)res.neutron_width>>31);
  res.B[1] = sqrt(abs(res.fission_width_1))*((int)res.fission_width_1>>31);
  res.B[2] = sqrt(abs(res.fission_width_2))*((int)res.fission_width_2>>31);
  //TODO:nomenclature
  res.ER = res.E - HALFI*res.radiation_width;
  res.EA = res.E - HALFI*res.total_width;
  res.COE1 = 0.0 - HALFI*res.neutron_width/res.sqrtE;
  res.COE2 = res.E - HALFI*res.absorption_width;
  penetration_shift_factor(iL, res);
    
}

void isotope::penetration_shift_factor(int iL, resonance &res){
  //rho rhohat treatment has been completed in initialize_l(),
  //channel_radius() now stores rho=kAP or ka
  double prho2 = pseudo_rho2[iL];
  double denominator;
  if(NOSHIFT){
    switch(iL){
    case 0:
      res.penetration_factor = 1.0;
      res.shift_factor = 0.0;
      res.COEF2 = 0.0 - res.COE1;
      break;
    case 1:
      res.penetration_factor = res.rho2/(1.0+res.rho2);
      res.shift_factor = 0.0;
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*prho2/res.penetration_factor;
    case 2:
      res.penetration_factor = 1.0/(9.0+3.0*res.rho2+res.rho4)*res.rho4;
      res.shift_factor = 0.0;
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*prho2*prho2/res.penetration_factor;
    case 3:
      res.penetration_factor = 1.0/(225.0+45.0*res.rho2+6.0*res.rho4 + res.rho6)*res.rho6;
      res.shift_factor = 0.0;
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*pow(prho2,3)/res.penetration_factor;
    default:
      break;
    }//end switch
  }
  else{
    switch(iL){
    case 0:
      res.penetration_factor = 1.0;
      res.shift_factor = 0.0;
      res.COEF2 = 0.0 - res.COE1;
      break;
    case 1:
      denominator = 1.0+res.rho2;
      res.penetration_factor = res.rho2/denominator;
      res.shift_factor = 0.0 - 1.0/denominator;
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*prho2/res.penetration_factor;
    case 2:
      denominator = 9.0+3.0*res.rho2+res.rho4;
      res.penetration_factor = 1.0/denominator*res.rho4;
      res.shift_factor = 0.0 - 1.0/denominator*(18.0+3.0*res.rho2);
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*prho2*prho2/res.penetration_factor;
    case 3:
      denominator = 225.0+45.0*res.rho2+6.0*res.rho4+res.rho6;
      res.penetration_factor = 1.0/denominator*res.rho6;
      res.shift_factor = 0.0 - 1.0/denominator*(675.0+90.8*res.rho2+6.0*res.rho4);
      res.COEF2 = 0.0 - 0.5*ONEI*res.neutron_width_0*pow(prho2,3)/res.penetration_factor;
    default:
      break;
    }//end switch

  }//end if else
  res.sqrtPF = sqrt(res.penetration_factor);
}

resonance &isotope::get_resonance(int iL, int iJ, int iR){
  return *(resonances[index(iL,iJ)]+iR);
}
int isotope::index(int iL, double j){
  int iJ = j - l_jdeg[2*iL+1]*0.5 ;
  int ind = iJ;
  for(int i=0;i<iL;i++)
    ind += l_jdeg[2*i];
  return ind;
}
int isotope::index(int iL, int iJ){
  int ind = iJ;
  for(int i=0;i<iL;i++)
    ind += l_jdeg[2*i];
  return ind;
}
void isotope::allocate_l(){
  l_jdeg = (unsigned*)malloc(number_l*sizeof(unsigned)*2);
  scattering_radius = (double*)malloc(number_l*sizeof(double));
  atomic_weight_ratio = (double*)malloc(number_l*sizeof(double));
  channel_radius = (double*)malloc(number_l*sizeof(double));
  pseudo_k = (double*)malloc(number_l*sizeof(double));
  pseudo_rho = (double*)malloc(number_l*sizeof(double));
  pseudo_rho2 = (double*)malloc(number_l*sizeof(double));
  factor = (double*)malloc(number_l*sizeof(double));
  pseudo_lambdabar2 = (double*)malloc(number_l*sizeof(double));
  pseudo_twolambdabar2 = (double*)malloc(number_l*sizeof(double));
}

void isotope::initialize_l(int iL){
  factor[iL] = atomic_weight_ratio[iL]/(atomic_weight_ratio[iL]+1.0);
  pseudo_lambdabar2[iL] = C2/factor[iL]/factor[iL];
  pseudo_twolambdabar2[iL] = 2.0*pseudo_lambdabar2[iL];
  pseudo_k[iL] = C1*factor[iL];
  if((0.0!=scattering_radius[iL])&&(1==flag_rcrs))
    channel_radius[iL] = scattering_radius[iL];
  else
    channel_radius[iL] = (1.23*pow(atomic_weight_ratio[iL],ONETRD) + 0.8)*0.1;
  pseudo_rho[iL]  = pseudo_k[iL]*channel_radius[iL];
  pseudo_rho2[iL] = pseudo_rho[iL] * pseudo_rho[iL];
}

int isotope::check_degeneracy(){
  //Note: all agular momentum numbers are doubled to integer here
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
  return sum;
}

void isotope::allocate_lj(int sum){
  number_resonances = (unsigned*)malloc(sum*sizeof(unsigned));
  number_channels   = (unsigned*)malloc(sum*sizeof(unsigned));
  channel_spin = (double*)malloc(sum*sizeof(double));
  gij          = (double*)malloc(sum*sizeof(double));
  resonances   = (resonance**)malloc(sum*sizeof(resonance*));
}

void isotope::initialize_lj(int iL, int iJ, int iLJ){
  channel_spin[iLJ] = 0.5*l_jdeg[2*iL+1] + iJ;
  gij[iLJ] = 0.5 * (2.0 * channel_spin[iLJ] + 1.0) / (2.0*target_spin + 1.0);
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

