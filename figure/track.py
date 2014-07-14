#!/bin/usr/python3
import numpy as np

data = np.loadtxt('../logfile/trackcnt')
colcnt = np.array([0,0])
wfucnt = np.array([0,0])
for i in range(len(data[:,0])):
  colcnt[(data[i][0]>2000.0)]+=1
  wfucnt[(data[i][0]>2000.0)]+=data[i][1]
print colcnt,wfucnt, wfucnt*1.0/colcnt
  
