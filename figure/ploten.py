#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

x=np.loadtxt('../boxtally',unpack=False);
plt.loglog(x[1000:len(x)-1]);
#plt.ylim([1e-5,1]);
plt.savefig('ASE.png');

