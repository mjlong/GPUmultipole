#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import ticker
import csv

from matplotlib import rc
rc('font',**{'family':'serif',
             'serif':['Palatino'],
             'size': 18})
rc('text',usetex=True)

def load_matrix_from_file(f):
    """
    This function is to load an ascii format matrix (float numbers separated by
    whitespace characters and newlines) into a numpy matrix object.
 
    f is a file object or a file path.
    """
 
    import types
    import numpy
 
    if type(f) == types.StringType:
        fo = open(f, 'r')
        matrix = load_matrix_from_file(fo)
        fo.close()
        return matrix
    elif type(f) == types.FileType:
        file_content = f.read().strip()
        file_content = file_content.replace('\r\n', ';')
        file_content = file_content.replace('\n', ';')
        file_content = file_content.replace('\r', ';')
 
        return numpy.matrix(file_content)
    
    raise TypeError('f must be a file object or a file name.')

mitd = load_matrix_from_file('../logfile/promitd')
energy = np.array(mitd[:,0])[:,0]
order = np.argsort(energy)
energy = energy[order]
sigt   = mitd[:,1][order]
siga   = mitd[:,2][order]
sigf   = mitd[:,3][order]
#fileList = ['proqgld','proqglf']#,'qtxf','qgld','qglf']
#titleist = ['QuickW --double','QuickW --float','QuickW texture --double', 'QuickW constant --double']

fileList = ['testfourierd','proqgld']
titleist = fileList;
i = 0
factor = math.sqrt(65536.0)
for stri in fileList:
    fileName = '../logfile/'+stri
    data = load_matrix_from_file('../logfile/'+stri)
    errt = abs(data[:,1][order] - sigt)/sigt*100
    tnorm = np.linalg.norm(abs(data[:,1][order] - sigt))
    erra = abs(data[:,2][order] - siga)/siga*100
    anorm = np.linalg.norm(abs(data[:,2][order] - siga))
    errf = abs(data[:,3][order] - sigf)/sigf*100
    fnorm = np.linalg.norm(abs(data[:,3][order] - sigf))
    plt.plot(energy,errt, energy,erra, energy,errf)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(('error sigt','error siga','error sigf'), loc=1,prop={'size':8})
    plt.title(titleist[i]+'compared with mitd')
    plt.xlabel('Energy(eV)')
    plt.ylabel('Relative error (\%)')
    plt.tight_layout()
    plt.savefig(titleist[i]+' compared with mitd_log.png')
    plt.figure()
    i = i + 1
    print tnorm, anorm, fnorm
    print stri+' done!\n'

