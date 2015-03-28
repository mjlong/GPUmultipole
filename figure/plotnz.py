#~/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import math

lines = np.loadtxt("t_nz",unpack=False);
plt.plot(lines[:,0],lines[:,1],lines[:,0],lines[:,2]);
plt.legend(('CPU','GPU'),loc=1,prop={'size':8})
plt.tight_layout();
plt.savefig('tvsnz.png')
