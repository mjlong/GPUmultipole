import numpy
def fac(n):
    return numpy.exp(numpy.sum(numpy.log(range(2,n+1))));
def lncnm(n,m):
    if(m==n):
        return 0;   
    if m<n/2:
        m = n-m;
    return numpy.sum(numpy.log(range(m+1,n+1)))-numpy.sum(numpy.log(range(2,n-m+1))); 

def cnm(n,m):
    if(m==n):
        return 1;   
    if m<n/2:
        m = n-m;
    return numpy.exp(numpy.sum(numpy.log(range(m+1,n+1)))-numpy.sum(numpy.log(range(2,n-m+1)))); 

def scnm(n,m):
    if(m==n):
        return 1;
    if(1==(n-m)):
        return n; 
    if(1==n):
        return 1;
    p=m*1.0/n;
    return numpy.sqrt(1.0/(m*(1-p)))/(p**m*(1-p)**(n-m));
