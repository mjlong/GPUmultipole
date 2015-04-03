#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import sys

x=np.loadtxt('../boxtally_'+sys.argv[1],delimiter=' ',unpack=False);
tally = x[:,0];
refer = x[:,1];

x=np.loadtxt('../ASE_EASE_'+sys.argv[1],delimiter=' ',unpack=False);
ASE = x[1:,0];
EASE= x[1:,1];
nhis = x[0,0];
nbat = x[0,1];
nbin = len(tally);

space = np.linspace(1,nbat,nbat);
plt.figure();
plt.plot(tally);
plt.plot(refer);
plt.savefig('dist'+'.png');

plt.figure();
plt.loglog(space,ASE,label='Average Square Error');
plt.loglog(space,EASE,label='Expectation of ASE');
plt.legend(loc='center left',prop={'size':10});
plt.xlabel('number of batch');
plt.ylabel('Absolute error (/cm/neutron)');
plt.title( 'numhist='+str(int(nhis)))
plt.savefig('ASE-n'+'.png');

x = np.loadtxt('../acc_0',delimiter=' ',unpack=False);
rho0 = x[0];
q    = x[1];
rhos = x[2:];
upto = len(x)-2;
plt.figure();
space = np.linspace(1,upto,upto);
plt.plot(space,rhos);
plt.plot(space,rho0*(q**space));
plt.savefig('acc.png');

