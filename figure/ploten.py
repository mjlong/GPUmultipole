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
plt.savefig(sys.argv[1]+'dist'+'.png');

plt.figure();
plt.loglog(space,ASE,label='Average Square Error: numhist='+str(nhis));
plt.loglog(space,EASE,label='Expectation of ASE: numhist='+str(nhis));
plt.legend(loc='lower left',prop={'size':16});
plt.xlabel('number of batch');
plt.ylabel('Absolute error (/cm/neutron)');
plt.savefig(sys.argv[1]+'ASE-n'+'.png');

