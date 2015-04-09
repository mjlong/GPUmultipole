# IPython log file

import matplotlib.pyplot as plt
import numpy
import h5py

files=h5py.File('Result_16384_4000_1000000strong.h5','r');
filew=h5py.File('../weak/Result_16384_1500_100000weak.h5','r');


batmeans=numpy.transpose(files['batchmeans']);
accmeans=numpy.transpose(files['batchaccumu']);
accmeanw=numpy.transpose(filew['batchaccumu']);

abat=accmeans.shape[1];
abatw=accmeanw.shape[1];

plt.plot(accmeans[:,-1]);
plt.plot(accmeanw[:,-1]);
#plt.ylim([0,0.5]);
plt.show()


width=files['width'][0][0]
meshes=files['num_cells'][0][0]
dx=width/meshes
batprob=accmeans*dx;
batprobw=accmeanw*dx;


for i in range(abat):
    batprob[:,i]=accmeans[:,i]/numpy.sum(accmeans[:,i]);

for i in range(abatw):
    batprobw[:,i]=accmeanw[:,i]/numpy.sum(accmeanw[:,i]);

# ============================================================================== 
# =====================N-th order Legendre moment===============================
# ============================================================================== 
space=numpy.linspace(-width/2+dx/2,width/2-dx/2,meshes)*2/width;
p1=space*1.0;
p2=1.0/2*(3.0*space**2-1.0);
p3=1.0/2*(5.0*space**3-3.0*space);
p4=1.0/8*(35.0*space**4-30.0*space**2+3);
plt.plot(p1);
plt.plot(p2);
plt.plot(p3);
plt.plot(p4);
plt.show();

xlog0=numpy.zeros([abat,1]);
xlog1=numpy.zeros([abat,1]);
xlog2=numpy.zeros([abat,1]);
xlog3=numpy.zeros([abat,1]);
xlog4=numpy.zeros([abat,1]);
for i in range(abat):
    xlog0[i] = 0.0-numpy.sum(    batprob[:,i]*numpy.log2(batprob[:,i]));
    xlog1[i] = 0.0-numpy.sum( p1*batprob[:,i]*numpy.log2(batprob[:,i]));
    xlog2[i] = 0.0-numpy.sum( p2*batprob[:,i]*numpy.log2(batprob[:,i]));
    xlog3[i] = 0.0-numpy.sum( p3*batprob[:,i]*numpy.log2(batprob[:,i]));
    xlog4[i] = 0.0-numpy.sum( p4*batprob[:,i]*numpy.log2(batprob[:,i]));
    
xlogw0=numpy.zeros([abatw,1]);
xlogw1=numpy.zeros([abatw,1]);
xlogw2=numpy.zeros([abatw,1]);
xlogw3=numpy.zeros([abatw,1]);
xlogw4=numpy.zeros([abatw,1]);
for i in range(abatw):
    xlogw0[i] = 0.0-numpy.sum(    batprobw[:,i]*numpy.log2(batprobw[:,i]));
    xlogw1[i] = 0.0-numpy.sum( p1*batprobw[:,i]*numpy.log2(batprobw[:,i]));
    xlogw2[i] = 0.0-numpy.sum( p2*batprobw[:,i]*numpy.log2(batprobw[:,i]));
    xlogw3[i] = 0.0-numpy.sum( p3*batprobw[:,i]*numpy.log2(batprobw[:,i]));
    xlogw4[i] = 0.0-numpy.sum( p4*batprobw[:,i]*numpy.log2(batprobw[:,i]));
 
# ==============================================================================   
# ==================== N-order entropy    ======================================
# ============================================================================== 
    
slogw2=numpy.zeros([abatw,1]);
slogw3=numpy.zeros([abatw,1]);
slogw4=numpy.zeros([abatw,1]);
for i in range(abatw):
    slogw2[i] =     numpy.sum( batprobwe[:,i]*(numpy.log2(batprobwe[:,i]))**2 );
    slogw3[i] = 0.0-numpy.sum( batprobwe[:,i]*(numpy.log2(batprobwe[:,i]))**3 );
    slogw4[i] =     numpy.sum( batprobwe[:,i]*(numpy.log2(batprobwe[:,i]))**4 );

slog2=numpy.zeros([abat,1]);
slog3=numpy.zeros([abat,1]);
slog4=numpy.zeros([abat,1]);
for i in range(abat):
    slog2[i] =     numpy.sum( batprob[:,i]*(numpy.log2(batprob[:,i]))**2 );
    slog3[i] = 0.0-numpy.sum( batprob[:,i]*(numpy.log2(batprob[:,i]))**3 );
    slog4[i] =     numpy.sum( batprob[:,i]*(numpy.log2(batprob[:,i]))**4 );
    
# ==============================================================================  
# ================ Plot  =======================================================        
# ==============================================================================        
plt.plot(slog4[0:9999]/(numpy.log2(meshes))**4);
plt.plot(numpy.ones([10000,1]));
plt.show();
plt.plot(slogw4[0:9999]/(numpy.log2(meshes))**4);
plt.plot(numpy.ones([10000,1]));
plt.show();
plt.plot(numpy.ones([10000,1]));
plt.plot(slogw2[2000:9999+2000]/(numpy.log2(meshes))**2);
plt.show();
plt.plot(numpy.ones([10000,1]));
plt.plot(slog2[2000:9999+2000]/(numpy.log2(meshes))**2);
plt.show();
plt.plot(numpy.ones([10000,1]));
plt.plot(slog2[2000:9999+2000]/(numpy.log2(meshes))**2);
ylim([0.8,1.002]);
plt.ylim([0.8,1.002]);
plt.show();
ylim([0.99,1.002]);
plt.ylim([0.99,1.002]);
plt.show();
plt.plot(slog2[2000:9999+2000]/(numpy.log2(meshes))**2);
plt.plot(slogw2[2000:9999+2000]/(numpy.log2(meshes))**2);
plt.plot(numpy.ones([10000,1]));
plt.show();
plt.plot(slog3[2000:9999+2000]/(numpy.log2(meshes))**3);
plt.plot(slogw3[2000:9999+2000]/(numpy.log2(meshes))**3);
plt.plot(numpy.ones([10000,1]));
plt.show();
plt.plot(slog4[2000:9999+2000]/(numpy.log2(meshes))**4);
plt.plot(slogw4[2000:9999+2000]/(numpy.log2(meshes))**4);
plt.plot(numpy.ones([10000,1]));
plt.show();
plt.plot(xlog4[2000:9999+2000]);
plt.plot(xlogw4[2000:9999+2000]);
plt.plot(numpy.zeros([10000,1]));
plt.show();
plt.plot(xlog4[2000:9999+2000],label='strong');
plt.plot(xlogw4[2000:9999+2000],label='weak');
plt.plot(numpy.zeros([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show();
plt.plot(xlog3[2000:9999+2000],label='strong');
plt.plot(xlogw3[2000:9999+2000],label='weak');
plt.plot(numpy.zeros([10000,1]),label='Reference');
plt.show();
plt.plot(xlog2[2000:9999+2000],label='strong');
plt.plot(xlogw2[2000:9999+2000],label='weak');
plt.plot(numpy.zeros([10000,1]),label='Reference');
plt.show();
plt.plot(xlog1[2000:9999+2000],label='strong');
plt.plot(xlogw1[2000:9999+2000],label='weak');
plt.plot(numpy.zeros([10000,1]),label='Reference');
plt.show();
plt.plot(xlog[2000:9999+2000],label='strong');
plt.plot(xlog2[2000:9999+2000],label='strong');
plt.plot(numpy.log2(meshes)*numpy.ones([10000,1]),label='Reference');
plt.show();
plt.plot(xlog[2000:9999+2000],label='strong');
plt.plot(xlogw[2000:9999+2000],label='weak');
plt.plot(numpy.log2(meshes)*numpy.ones([10000,1]),label='Reference');
plt.show();
plt.plot(xlog[2000:9999+2000],label='strong');
plt.plot(xlogw[2000:9999+2000],label='weak');
plt.plot(numpy.log2(meshes)*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show();
plt.plot(xlog1[2000:9999+2000],label='strong');
plt.plot(xlogw1[2000:9999+2000],label='weak');
plt.plot(numpy.log2(meshes)*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show();
plt.plot(xlog1[2000:9999+2000],label='strong');
plt.plot(xlogw1[2000:9999+2000],label='weak');
plt.plot(0*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show();
plt.plot(xlog1[2000:9999+2000],label='strong');
plt.plot(xlogw1[2000:9999+2000],label='weak');
plt.plot(0*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='upper right');
plt.show();
plt.plot(xlog2[2000:9999+2000],label='strong');
plt.plot(xlogw2[2000:9999+2000],label='weak');
plt.plot(0*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='upper right');
plt.show();
plt.plot(xlog3[2000:9999+2000],label='strong');
plt.plot(xlogw3[2000:9999+2000],label='weak');
plt.plot(0*numpy.ones([10000,1]),label='Reference');
plt.legend(loc='upper right');
plt.show();
plt.plot(slog3[0:20000]);
plt.plot(slogw3[0:20000]);
plt.show();
plt.plot(slog2[2000:9999+2000]/numpy.log2(meshes)**2,label='strong');
plt.plot(slogw2[2000:9999+2000]/numpy.log2(meshes)**2,label='weak');
plt.plot(numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show()
plt.plot(slog2[0:20000]);
plt.plot(slog3[0:20000]);
plt.show()
plt.plot(slog2[0:20000]/numpy.log2(meshes)**2);
plt.plot(slog3[0:20000]/numpy.log2(meshes)**3);
plt.plot(slog4[0:20000]/numpy.log2(meshes)**4);
plt.show()
plt.plot(slog3[2000:9999+2000]/numpy.log2(meshes)**3,label='strong');
plt.plot(slogw3[2000:9999+2000]/numpy.log2(meshes)**3,label='weak');
plt.plot(numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show()
plt.plot(slog4[2000:9999+2000]/numpy.log2(meshes)**4,label='strong');
plt.plot(slogw4[2000:9999+2000]/numpy.log2(meshes)**4,label='weak');
plt.plot(numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show()
plt.plot(xlog[2000:9999+2000]/numpy.log2(meshes),label='strong');
plt.plot(xlogw[2000:9999+2000]/numpy.log2(meshes),label='weak');
plt.show()
plt.plot(xlog[2000:9999+2000]/numpy.log2(meshes),label='strong');
plt.plot(xlogw[2000:9999+2000]/numpy.log2(meshes),label='weak');
plt.plot(numpy.ones([10000,1]),label='Reference');
plt.legend(loc='lower right');
plt.show()
