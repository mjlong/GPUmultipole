#include "process.h"
void cnt2flux(MemStruct HostMem, unsigned numhis, float dx, unsigned meshes, unsigned nbat){
  int im, ib;
  for(ib=0;ib<nbat;ib++){
    for(im=0;im<meshes;im++){
      HostMem.accmeans[ib*meshes+im] = HostMem.acccnt[ib*meshes+im]/(numhis*dx*(ib+1));
    }
  }

  ib = 0;
  for(im=0;im<meshes;im++)
    HostMem.batchmeans[im] = HostMem.accmeans[im];

  for(ib=1;ib<nbat;ib++){
    for(im=0;im<meshes;im++){
      HostMem.batchmeans[ib*meshes+im] = (HostMem.acccnt[ib*meshes+im]-HostMem.acccnt[(ib-1)*meshes+im])/(numhis*dx);
    }
  }

  
}

void getASE(double *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, double ref, double* ASE){
  int ib,im,index;
  for(ib=ubat;ib<nbat;ib++)
    ASE[ib-ubat]=0.0;
  for(ib=ubat;ib<nbat;ib++){
    for(im=0;im<meshes;im++){
      index = ib*meshes+im;
      ASE[ib-ubat] += (accmeans[index]-ref)*(accmeans[index]-ref)/meshes;
    }
  }
}


void getCOR(double *batmeans, unsigned meshes, unsigned nbat, unsigned ubat, unsigned upto,double *COR){
  int im,ib;
  double *batmeans_active = batmeans+ubat*meshes;
  for(im=0;im<meshes;im++){
    for(ib=0;ib<upto;ib++){
      COR[im*upto+ib] = autok(batmeans_active,nbat-ubat,ib+1,meshes,im);
    }
  }
  

}

double autok(double *batmeans, unsigned n, unsigned k, unsigned meshes, unsigned im){
  double sum1=0; double sum2=0; double sum3=0; double sum4=0; double sum5 = 0;
  double xi,xik;
  int ib;
  for(ib=0;ib<n-k;ib++){
    //printf("in autok, ib=%d\n",ib);
    xi = batmeans[ib*meshes+im];
    xik= batmeans[(ib+k)*meshes+im];
    sum1+=xi;
    sum2+=xik;
    sum3+=xi*xik;
    sum4+=xi*xi;
    sum5+=xik*xik;
  }
  
  return ((n-k)*sum3-sum1*sum2)/sqrt(((n-k)*sum4-sum1*sum1)*((n-k)*sum5-sum2*sum2));
}

void fitrho(double* rho, unsigned m, double* rho0, double* q){
  double sum1=0; double sum2=0;
  for(int i=0;i<m;i++){
    sum1 += log(rho[i]);
    sum2 += log(rho[i])*(i+1);
  }
  *rho0 = exp(2.0*(sum1+2.0*m*sum1-3.0*sum2)/(m*(m-1.0)));
  *q = exp((6.0*sum1+6.0*m*sum1-12.0*sum2)/(m*(1.0-m*m)));
  /*
    m = length(ρ);
    A = [m              0.5*m*(m+1);
         0.5*m*(m+1)    1.0/6*m*(m+1)*(2*m+1)];
    c = inv(A)*[sum(log(ρ)), sum(log(ρ).*range(1,m))];
    return exp(c[1]),exp(c[2]); #ρ0, q
  */
}

void fitrho1(double* rho, unsigned m, double* rho0, double* q){
  double sum1 = 0;
  double rho1 = rho[0];
  for(int i=0;i<m;i++)
    sum1 += log(rho[i])*i;
  *q = exp((sum1-m*(m-1.0)*0.5*log(rho1))/(1.0/6*(m-1.0)*m*(2.0*m-1.0)));
  *rho0 = rho1/(*q);
}


void fitall(double *rhos,unsigned upto, unsigned meshes, double *rho0s, double *qs){
  for(int im=0;im<meshes;im++){
    fitrho(rhos+im*upto, upto, rho0s+im, qs+im);
  }
}

void fitall1(double *rhos,unsigned upto, unsigned meshes, double *rho0s, double *qs){
  for(int im=0;im<meshes;im++){
    fitrho1(rhos+im*upto, upto, rho0s+im, qs+im);
  }
}
