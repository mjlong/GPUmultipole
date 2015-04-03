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

void getASE(float *accmeans,unsigned meshes, unsigned nbat, unsigned ubat, float ref, float* ASE){
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
