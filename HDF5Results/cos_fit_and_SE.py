# IPython log file

import numpy
import h5py
import matplotlib.pyplot as plt
pi=numpy.pi;
sin = numpy.sin;
cos = numpy.cos;
exp = numpy.exp;
log = numpy.log;
log2= numpy.log2;


files=h5py.File('cosstrong/Result_16384_4000_1000000.h5','r');
filew=h5py.File('cosweak/Result_16384_4000_1000000.h5','r');
meshes=files['num_cells'][0][0]
width=files['width'][0][0]
dx=width/meshes;
accmeans=numpy.transpose(files['batchaccumu']);
accmeanw=numpy.transpose(filew['batchaccumu']);
plt.plot(accmeans[:,-1])
plt.plot(accmeanw[:,-1])
plt.show()
batmeans=numpy.transpose(files['batchmeans']);
batmeanw=numpy.transpose(filew['batchmeans']);
plt.plot(batmeans[:,-1])
plt.plot(accmeans[:,-1])
plt.show()
plt.plot(batmeanw[:,-1])
plt.plot(accmeanw[:,-1])
plt.show()
batprobs=accmeans*dx;
batprobw=accmeanw*dx;
abat = accmeans.shape[1]
for i in range(abat):
    batprobs[:,i]=accmeans[:,i]/numpy.sum(accmeans[:,i]);
    batprobw[:,i]=accmeanw[:,i]/numpy.sum(accmeanw[:,i]);
    
batprobdens=batprobs/dx;
batprobdenw=batprobw/dx;
plt.show()
xspace=numpy.linspace(-width/2+dx/2,width/2-dx/2,meshes);
plt.plot(xspace,batprobdens[:,-1])
plt.show();

As=numpy.max(batprobdens[:,-1]);
Ls=numpy.pi*xspace[0]/numpy.arccos(0.5*(batprobdens[0,-1]+batprobdens[-1,-1])/As)
fits=As*numpy.cos(numpy.pi*xspace/Ls);
plt.plot(xspace,fits);
plt.plot(xspace,batprobdens[:,-1]);
plt.show();

Aw=numpy.max(batprobdenw[:,-1]);
Lw=numpy.pi*xspace[0]/numpy.arccos(0.5*(batprobdenw[0,-1]+batprobdenw[-1,-1])/Aw);
fitw=Aw*numpy.cos(numpy.pi*xspace/Lw);
plt.plot(xspace,batprobdenw[:,-1]);
plt.plot(xspace,fitw);
plt.show();


SEs = calculateSE(batprobs);
plt.plot(SEs);
plt.show();

plt.plot(SEs[0:9999]);
plt.plot(numpy.ones([10000,1])*numpy.log2(meshes));
plt.show();

tSEs = SEcos(Ls,dx,width);
plt.plot(SEs[0:9999]);
plt.plot(numpy.ones([10000,1])*tSEs);
plt.show();

fitPs=numpy.cos(pi*xspace/Ls)*numpy.sin(dx*pi*0.5/Ls)/numpy.sin(width*pi*0.5/Ls);
tSEs_cal = numpy.sum(fitPs*numpy.log2(fitPs))*(-1);
plt.plot(SEs[0:9999]);
plt.plot(numpy.ones([10000,1])*tSEs,label='theoretical SE');
plt.plot(numpy.ones([10000,1])*tSEs_cal,label='theoretical SE (sum)');
plt.legend(loc='lower center');
plt.show();
plt.plot(SEs[0:9999],label='SE');
plt.plot(numpy.ones([10000,1])*tSEs_cal,label='theoretical SE (sum)');
plt.plot(numpy.ones([10000,1])*tSEs,label='theoretical SE (analytic integration)');
plt.legend(loc='lower center');
plt.show();
tSEw=SEcos(Lw,dx,width);
fitPw=numpy.cos(pi*xspace/Lw)*numpy.sin(dx*pi*0.5/Lw)/numpy.sin(width*pi*0.5/Lw);
tSEw_cal = numpy.sum(fitPw*numpy.log2(fitPw))*(-1);
SEw = calculateSE(batprobw);
plt.plot(SEw[0:9999],label='SE');
plt.plot(numpy.ones([10000,1])*tSEw,label='theoretical SE (analytic integration)');
plt.plot(numpy.ones([10000,1])*tSEw_cal,label='theoretical SE (sum)');
plt.legend(loc='lower center');
plt.show();
plt.plot(SEs[0:9999],label='SE');
plt.plot(numpy.ones([10000,1])*tSEs,label='theoretical SE (analytic integration)');
plt.plot(numpy.ones([10000,1])*tSEs_cal,label='theoretical SE (sum)');
plt.legend(loc='lower center');
plt.show();

nhis=files['num_history'][0][0]
nneu = int(numpy.sum(accmeans[:,-1]*dx*nhis))
nneu = nhis;

Nb=8;
nn=100;
Ltemp=3.0;
atemp=2.0;
htemp=atemp/Nb;
xtemp=numpy.linspace(-atemp*0.5+htemp*0.5,atemp*0.5-htemp*0.5,Nb);
xtemp
ptemp=numpy.cos(pi*xtemp/Ltemp)*numpy.sin(pi*htemp*0.5/Ltemp)/numpy.sin(atemp*pi*0.5/Ltemp);
ESEcos(Nb,nn,xtemp,ptemp)

