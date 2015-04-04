#!/usr/bin/python3
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

num = len(sys.argv)-1;#num must be at least 2
print num;
dir = '../rand1/'
colors = ('b','g','r','c','m','y');
#'''
# ==============================================================================
# =============plot ASME together, use EASME from largest batch ================
print 'plotting ASME together'
plt.figure();
for i in range(num):
  print 'plotting ASME of'+sys.argv[i+1]
  x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[i+1],delimiter=' ',unpack=False);
  ASE = x[1:,0];
  #EASE= x[1:,1];
  nhis = x[0,0];
  nbat = x[0,1];
  print nhis, nbat
  #nbin = len(tally);
  space = np.linspace(1,nbat,nbat);
  plt.loglog(space,ASE,color=colors[i],label='Average Square Mean Error(ASME) nhis:'+str(int(nhis)));
  
print 'adding EASE from'+sys.argv[num];
x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[num],delimiter=' ',unpack=False);
EASE= x[1:,1];

plt.loglog(space,EASE,'-',color='black',label='Expectation of ASME');
plt.loglog(space,EASE[0]/space,'-.',color='black', label='Variance over N');
plt.legend(loc='lower left',prop={'size':10});
plt.xlabel('number of batch');
plt.ylabel('Absolute error (/cm/neutron)^2');
plt.savefig('ASME-n'+'.png');
#'''
# ==============================================================================
# =======================plotting all batches in one============================
print 'plotting all batches in one'
plt.figure();
for i in range(num):
  x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[i+1],delimiter=' ',unpack=False);
  ASE = x[1:,0];
  EASE= x[1:,1];
  nhis = x[0,0];
  nbat = x[0,1];
  #nbin = len(tally);
  space = np.linspace(1,nbat,nbat);
  plt.loglog(space,ASE,'.',color=colors[i],label='Average Square Mean Error(ASME) nhis:'+str(int(nhis)));
  plt.loglog(space,EASE/EASE[0]*ASE[0],color=colors[i],label='shifted Expectation of ASME:'+str(int(nhis)));
  print 'plotting ASME of'+sys.argv[i+1]

plt.legend(loc='lower left',prop={'size':10});
plt.xlabel('number of histories');
plt.ylabel('Absolute error (/cm/neutron)^2');
plt.savefig('ASME-n-more'+'.png');

# ==============================================================================
# =======================plotting batches separately============================
#'''
for i in range(num):
  print 'plotting separate ASE EASE '+sys.argv[i+1]
  x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[i+1],delimiter=' ',unpack=False);
  ASE = x[1:,0];
  EASE= x[1:,1];
  nhis = x[0,0];
  nbat = x[0,1];
  plt.figure();
  space = np.linspace(1,nbat,nbat);
  plt.loglog(space,ASE,label='Average Square Mean Error (ASME)');
  plt.loglog(space,EASE,'-',color='gray',label='Expectation of ASME (EASME)');
  plt.loglog(space,EASE/EASE[0]*ASE[0],'-.',color='black',label='shifted EASME');
  plt.loglog(space,EASE[0]/space,'-',color='gray', label='Variance over N');
  plt.loglog(space,ASE[0]/space,'-.',color='black', label='shifted Variance over N');

  plt.legend(loc='lower left',prop={'size':8});
  plt.xlabel('number of batch');
  plt.ylabel('Absolute error (/cm/neutron)^2');
  plt.title( 'numhist='+str(int(nhis)))
  plt.savefig('ASME-n'+sys.argv[i+1]+'.png');
#'''
'''
  plt.figure();
  print 'plotting independent distribution from'+sys.argv[i+1]
  x=np.loadtxt(dir+'boxtally_'+sys.argv[i+1],delimiter=' ',unpack=False);
  tally = x[:,0];
  refer = x[:,1];
  plt.figure();
  plt.plot(tally);
  plt.plot(refer);
  plt.savefig('dist'+sys.argv[i+1]+'.png');
'''
  
'''
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

'''
