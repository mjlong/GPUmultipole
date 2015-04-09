def calculateSE(batprob):
    abat=batprob.shape[1];
    SEs=numpy.zeros([abat,1]);
    for i in range(abat):
        SEs[i] = numpy.sum(batprob[:,i]*numpy.log2(batprob[:,i]));
    return SEs*(-1);

def SEcos(L,h,a):
    x=numpy.pi*a/(4*L);
    sin4x=numpy.sin(x);
    cos4x=numpy.cos(x);
    cos2x=numpy.cos(2*x);
    sin2x=numpy.sin(2*x);
    return -numpy.log2(numpy.pi*h/(2*L*sin2x))-(numpy.log2((cos4x+sin4x)/(cos4x-sin4x))+sin2x*(-numpy.log2(numpy.e)+numpy.log2(cos2x)))/sin2x;

 
import fastfactorial
import numpy 
def ESEany(Nb,nn,space,ps):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    for ii in range(Nb):
        sum = sum + exp(                                                  +log(ps[ii])*nn +log(nn)   + log(log2(nn)));#xi=nn
        sum = sum + exp(                 log(nn)      + log(1-ps[ii]) +log(ps[ii])*(nn-1) +log(nn-1) + log(log2(nn-1)));#xi=nn-1
        for xi in range(2,(nn+1)/2):
            lncnx = fastfactorial.lncnm(nn,xi);
            sum = sum + exp(lncnx+log(1-ps[ii])*(nn-xi)+log(ps[ii])*xi+log(xi)+log(log2(xi)))+exp(lncnx+log(1-ps[ii])*xi +log(ps[ii])*(nn-xi)+log(nn-xi)+log(log2(nn-xi)));
        if(0==nn%2): 
            xi = nn/2;
            sum = sum + exp(fastfactorial.lncnm(nn,xi)+ log(1-ps[ii])*(nn-xi) +log(ps[ii])*xi +log(xi) + log(log2(xi)));
    return log2(nn)-sum/nn;
    

def ESEany(Nb,nn,space,ps,nm):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    for ii in range(Nb):
        sum = sum + exp(                                                  +log(ps[ii])*nn +log(nn)   + log(log2(nn)));#xi=nn
        sum = sum + exp(                 log(nn)      + log(1-ps[ii]) +log(ps[ii])*(nn-1) +log(nn-1) + log(log2(nn-1)));#xi=nn-1
        for xi in range(2,(nn+1)/2):
            lncnx = fastfactorial.lncnm(nn,xi);
            sum = sum + exp(lncnx+log(1-ps[ii])*(nn-xi)+log(ps[ii])*xi+log(xi)+log(log2(xi)))+exp(lncnx+log(1-ps[ii])*xi +log(ps[ii])*(nn-xi)+log(nn-xi)+log(log2(nn-xi)));
        if(0==nn%2): 
            xi = nn/2;
            sum = sum + exp(fastfactorial.lncnm(nn,xi)+ log(1-ps[ii])*(nn-xi) +log(ps[ii])*xi +log(xi) + log(log2(xi)));
    return log2(nn)-sum/nn;

def BnpMoment(n,p,m):
    if(2==m):
    if(3==m):
    if(4==m):
