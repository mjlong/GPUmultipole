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
