#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys

num = len(sys.argv)-1;#num must be at least 2
print num;
dir = '../rand1/'
colors = ('b','g','r','c','m','y');

# ==============================================================================
# ============================plot over batches ================================
print 'plotting over batches'
plt.figure();
for i in range(num):
  print 'plotting ASE of'+sys.argv[i+1]
  x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[i+1],delimiter=' ',unpack=False);
  ASE = x[1:,0];
  #EASE= x[1:,1];
  nhis = x[0,0];
  nbat = x[0,1];
  print nhis, nbat
  #nbin = len(tally);
  space = np.linspace(1,nbat,nbat);
  plt.loglog(space,ASE,color=colors[i],label='Average Square Error nhis:'+str(int(nhis)));
  
print 'plotting EASE from'+sys.argv[num];
x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[num],delimiter=' ',unpack=False);
EASE= x[1:,1];

plt.loglog(space,EASE,'-',color='black',label='Expectation of ASE');
plt.loglog(space,EASE[0]/space,'-.',color='black', label='Variance over N');
plt.legend(loc='lower left',prop={'size':10});
plt.xlabel('number of batch');
plt.ylabel('Absolute error (/cm/neutron)^2');
plt.savefig('ASE-n'+'.png');

# ==============================================================================
# ============================plot over histories ================================
#'''
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
  plt.loglog(space,ASE,'.',color=colors[i],label='Average Square Error nhis:'+str(int(nhis)));
  plt.loglog(space,EASE,color=colors[i],label='Expectation of ASE:'+str(int(nhis)));
  print 'plotting ASE of'+sys.argv[i+1]

plt.legend(loc='lower left',prop={'size':10});
plt.xlabel('number of histories');
plt.ylabel('Absolute error (/cm/neutron)^2');
plt.savefig('ASE-n-more'+'.png');
#'''

'''
for i in range(num):
  print 'plotting independent ASE EASE '+sys.argv[i+1]
  x=np.loadtxt(dir+'ASE_EASE_'+sys.argv[i+1],delimiter=' ',unpack=False);
  ASE = x[1:,0];
  EASE= x[1:,1];
  nhis = x[0,0];
  nbat = x[0,1];
  plt.figure();
  space = np.linspace(1,nbat,nbat);
  plt.loglog(space,ASE,label='Average Square Error');
  plt.loglog(space,EASE,color='black',label='Expectation of ASE');
  plt.loglog(space,EASE[0]/space,'-.',color='black', label='Variance over N');
  plt.legend(loc='lower left',prop={'size':10});
  plt.xlabel('number of batch');
  plt.ylabel('Absolute error (/cm/neutron)^2');
  plt.title( 'numhist='+str(int(nhis)))
  plt.savefig('ASE-n'+sys.argv[i+1]+'.png');

  plt.figure();
  print 'plotting independent distribution from'+sys.argv[i+1]
  x=np.loadtxt(dir+'boxtally_'+sys.argv[i+1],delimiter=' ',unpack=False);
  tally = x[:,0];
  refer = x[:,1];
  plt.figure();
  plt.plot(tally);
  plt.plot(refer);
  plt.savefig('dist'+sys.argv[i+1]+'.png');
  

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
