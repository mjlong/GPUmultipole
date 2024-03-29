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

def ESEunilec(Nb,nn):
    return numpy.log2(Nb)+numpy.log2(1-(1-1.0/Nb)**nn);
def ESEuni1(Nb,nn,ps):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    sum = sum + exp(                                                  +log(ps)*nn +log(nn)   + log(log2(nn)));#xi=nn
    sum = sum + exp(                 log(nn)      + log(1-ps) +log(ps)*(nn-1) +log(nn-1) + log(log2(nn-1)));#xi=nn-1
    for xi in range(2,(nn+1)/2):
        lncnx = fastfactorial.lncnm(nn,xi);
        sum = sum + exp(lncnx+log(1-ps)*(nn-xi)+log(ps)*xi+log(xi)+log(log2(xi)))+exp(lncnx+log(1-ps)*xi +log(ps)*(nn-xi)+log(nn-xi)+log(log2(nn-xi)));
    if(0==nn%2): 
        xi = nn/2;
        sum = sum + exp(fastfactorial.lncnm(nn,xi)+ log(1-ps)*(nn-xi) +log(ps)*xi +log(xi) + log(log2(xi)));
    return log2(nn)-Nb*sum/nn;

def ESEany1(Nb,nn,ps):
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
    

def ESEuni2(Nb,nn,ps,eps):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    np = int(nn*ps);
    n1 = max(2,np-int(nn**(0.5+eps)));
    n2 = min(np+int(nn**(0.5+eps))+1,nn+1);
    sum = 0.0;
    for xi in range(n1,n2):
        lncnx = fastfactorial.lncnm(nn,xi);
        sum = sum + exp(lncnx+log(1-ps)*(nn-xi)+log(ps)*xi+log(xi)+log(log2(xi)));
    return log2(nn)-Nb*sum/nn;

def ESEany2(Nb,nn,ps,eps):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    for ii in range(Nb):
        np = int(nn*ps[ii]);
        n1 = max(2,np-int(nn**(0.5+eps)));
        n2 = min(np+int(nn**(0.5+eps))+1,nn+1);
        for xi in range(n1,n2):
            lncnx = fastfactorial.lncnm(nn,xi);
            sum = sum + exp(lncnx+log(1-ps[ii])*(nn-xi)+log(ps[ii])*xi+log(xi)+log(log2(xi)));
    return log2(nn)-sum/nn;
 



def ESEany(Nb,nn,ps,nm):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    for ii in range(Nb):
        np = nn*ps[ii];
        sum = sum + np*log(np);
        for im in range(2,nm+1):
            sum = sum + BnpMoment(nn,ps[ii],im)/(im*(im-1)*np**(im-1))*(-1)**im;    
    return log2(nn)-sum/log(2.0)/nn;

def ESEuni3(Nb,nn,ps):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    np = nn*ps;
    sum = sum + np*log(np)+0.5*(1-ps);
    return log2(nn)-sum*Nb/log(2.0)/nn;


def ESEuni(Nb,nn,ps,nm):
    log = numpy.log
    log2= numpy.log2
    exp = numpy.exp
    sum = 0.0;
    np = nn*ps;
    sum = sum + np*log(np);
    for im in range(2,nm+1):
        sum = sum + BnpMoment(nn,ps,im)/(im*(im-1)*np**(im-1))*(-1)**im;    
    return log2(nn)-sum*Nb/log(2.0)/nn;

def BnpMoment(n,P,m):
    if(2==m):
        return -n*(-1+P)*P
    if(3==m):
        return n*(-1+P)*P*(-1+2*P)
    if(4==m):
        return n*(-1+P)*P*(-1+3*(-2+n)*(-1+P)*P)
    if(5==m):
        return -n*(-1+P)*P*(-1+2*P)*(-1+2*(-6+5*n)*(-1+P)*P)
    if(6==m):
        return -n*(-1+P)*P*(1+5*(-1+P)*P*(6*(1-2*P)**2+3*n**2*(-1+P)*P+n*(-5-26*(-1+P)*P)))
    if(7==m):
        return n*(-1+P)*P*(-1+2*P)*(1+(-1+P)*P*(105*n**2*(-1+P)*P+60*(1+6*(-1+P)*P)-14*n*(4+33*(-1+P)*P)))
    if(8==m):
        return n*(-1+P)*P*(-1+7*(-1+P)*P*(15*n**3*(-1+P)**2*P**2-10*n**2*(-1+P)*P*(7+34*(-1+P)*P)-6*(3+40*(-1+P)*P*(1+3*(-1+P)*P))+n*(17+4*(-1+P)*P*(77+261*(-1+P)*P))))
    if(9==m):
        return -n*(-1+P)*P*(-1+2*P)*(-1+2*(-1+P)*P*(630*n**3*(-1+P)**2*P**2-7*n**2*(-1+P)*P*(137+944*(-1+P)*P)-126*(1+20*(1-2*P)**2*(-1+P)*P)+3*n*(41+6*(-1+P)*P*(193+892*(-1+P)*P))))
    if(10==m):
        return -n*(-1+P)*P*(1+3*(-1+P)*P*(315*n**4*(-1+P)**3*P**3-1050*n**3*(-1+P)**2*P**2*(3+14*(-1+P)*P)+35*n**2*(-1+P)*P*(65+4*(-1+P)*P*(236+723*(-1+P)*P))+10*(17+84*(-1+P)*P*(7+6*P*(-3+2*P))*(1+6*P*(-1+2*P)))+n*(-167-2*(-1+P)*P*(4073+36*(-1+P)*P*(1115+2886*(-1+P)*P)))))
    if(11==m):
        return n*(-1+P)*P*(-1+2*P)*(1+(-1+P)*P*(17325*n**4*(-1+P)**3*P**3-1540*n**3*(-1+P)**2*P**2*(37+229*(-1+P)*P)+165*n**2*(-1+P)*P*(139+4*(-1+P)*P*(706+2787*(-1+P)*P))+60*(17+126*(-1+P)*P*(7+80*(-1+P)*P*(1+3*(-1+P)*P)))-22*n*(46+9*(-1+P)*P*(383+40*(-1+P)*P*(128+419*(-1+P)*P)))))
    if(12==m):
        return n*(-1+P)*P*(-1+11*(-1+P)*P*(945*n**5*(-1+P)**4*P**4-1575*n**4*(-1+P)**3*P**3*(11+50*(-1+P)*P)+35*n**3*(-1+P)**2*P**2*(787+4*(-1+P)*P*(2519+7207*(-1+P)*P))-12*n**2*(-1+P)*P*(563+(-1+P)*P*(18589+3*(-1+P)*P*(51095+119506*(-1+P)*P)))-6*(31+40*(-1+P)*P*(64+63*(-1+P)*P*(19+120*(-1+P)*P*(1+2*(-1+P)*P))))+n*(185+4*(-1+P)*P*(5528+9*(-1+P)*P*(13411+20*(-1+P)*P*(4609+9722*(-1+P)*P))))))
    if(13==m):
        return -n*(-1+P)*P*(-1+2*P)*(-1+(-1+P)*P*(270270*n**5*(-1+P)**4*P**4-15015*n**4*(-1+P)**3*P**3*(109+628*(-1+P)*P)+2860*n**3*(-1+P)**2*P**2*(520+(-1+P)*P*(8737+31002*(-1+P)*P))-1716*n**2*(-1+P)*P*(137+(-1+P)*P*(6250+3*(-1+P)*P*(22015+62892*(-1+P)*P)))-132*(31+240*(-1+P)*P*(16+21*(-1+P)*P*(19+30*(-1+P)*P*(5+12*(-1+P)*P))))+26*n*(157+24*(-1+P)*P*(1189+3*(-1+P)*P*(11687+5*(-1+P)*P*(20275+51684*(-1+P)*P))))))
    if(14==m):
        return -n*(-1+P)*P*(1+13*(-1+P)*P*(10395*n**6*(-1+P)**5*P**5-24255*n**5*(-1+P)**4*P**4*(13+58*(-1+P)*P)+5390*n**4*(-1+P)**3*P**3*(173+(-1+P)*P*(2041+5558*(-1+P)*P))-308*n**3*(-1+P)**2*P**2*(1727+(-1+P)*P*(45856+(-1+P)*P*(335855+729402*(-1+P)*P)))+210*(3+44*(-1+P)*P*(1+12*(-1+P)*P)*(13+108*(-1+P)*P*(3+20*(-1+P)*P*(1+2*(-1+P)*P))))+77*n**2*(-1+P)*P*(731+8*(-1+P)*P*(6163+9*(-1+P)*P*(11499+(-1+P)*P*(68133+130688*(-1+P)*P))))+n*(-629-2*(-1+P)*P*(88201+144*(-1+P)*P*(26735+(-1+P)*P*(369593+5*(-1+P)*P*(391937+699612*(-1+P)*P)))))))
    if(15==m):
        return n*(-1+P)*P*(-1+2*P)*(1+(-1+P)*P*(4729725*n**6*(-1+P)**5*P**5-210210*n**5*(-1+P)**4*P**4*(226+1237*(-1+P)*P)+50050*n**4*(-1+P)**3*P**3*(1627+2*(-1+P)*P*(12067+39635*(-1+P)*P))-40040*n**3*(-1+P)**2*P**2*(773+(-1+P)*P*(26678+3*(-1+P)*P*(80380+207929*(-1+P)*P)))+91*n**2*(-1+P)*P*(24751+8*(-1+P)*P*(285440+9*(-1+P)*P*(679193+2*(-1+P)*P*(2451139+5548328*(-1+P)*P))))+16380*(1+22*(-1+P)*P*(13+160*(-1+P)*P*(4+63*(-1+P)*P*(1+6*(-1+P)*P*(1+2*(-1+P)*P)))))-6*n*(2728+3*(-1+P)*P*(385387+80*(-1+P)*P*(282972+7*(-1+P)*P*(703943+6*(-1+P)*P*(750983+1571266*(-1+P)*P)))))))
    if(16==m):
        return n*(-1+P)*P*(-1+(-1+P)*P*(2027025*n**7*(-1+P)**6*P**6-18918900*n**6*(-1+P)**5*P**5*(5+22*(-1+P)*P)+210210*n**5*(-1+P)**4*P**4*(2221+20*(-1+P)*P*(1237+3245*(-1+P)*P))-10010*n**4*(-1+P)**3*P**3*(51037+2*(-1+P)*P*(589063+4*(-1+P)*P*(990875+2033364*(-1+P)*P)))+5005*n**3*(-1+P)**2*P**2*(26929+8*(-1+P)*P*(165793+(-1+P)*P*(2344229+6*(-1+P)*P*(2079290+3707031*(-1+P)*P))))-78*n**2*(-1+P)*P*(88201+2*(-1+P)*P*(5750849+6*(-1+P)*P*(29372015+2*(-1+P)*P*(166333489+40*(-1+P)*P*(19419148+31752267*(-1+P)*P)))))-6*(5461+10920*(-1+P)*P*(217+33*(-1+P)*P*(457+240*(-1+P)*P*(43+14*(-1+P)*P*(29+60*(-1+P)*P*(2+3*(-1+P)*P))))))+3*n*(10917+4*(-1+P)*P*(1758119+3*(-1+P)*P*(48623879+40*(-1+P)*P*(30332101+14*(-1+P)*P*(21801607+60*(-1+P)*P*(1571266+2434119*(-1+P)*P))))))))
    if(17==m):
        return -n*(-1+P)*P*(-1+2*P)*(-1+2*(-1+P)*P*(45945900*n**7*(-1+P)**6*P**6-3573570*n**6*(-1+P)**5*P**5*(199+1048*(-1+P)*P)+170170*n**5*(-1+P)**4*P**4*(12059+50*(-1+P)*P*(3275+10144*(-1+P)*P))-85085*n**4*(-1+P)**3*P**3*(17921+4*(-1+P)*P*(129061+2*(-1+P)*P*(521485+1250424*(-1+P)*P)))+3094*n**3*(-1+P)**2*P**2*(92947+2*(-1+P)*P*(2950595+6*(-1+P)*P*(8543131+2*(-1+P)*P*(27005333+55806276*(-1+P)*P))))-51*n**2*(-1+P)*P*(204911+4*(-1+P)*P*(9067454+3*(-1+P)*P*(117631295+4*(-1+P)*P*(404777303+80*(-1+P)*P*(27842759+52444242*(-1+P)*P)))))-6*(5461+16380*(-1+P)*P*(217+44*(-1+P)*P*(457+60*(-1+P)*P*(215+84*(-1+P)*P*(29+20*(-1+P)*P*(7+12*(-1+P)*P))))))+17*n*(1927+18*(-1+P)*P*(103847+4*(-1+P)*P*(2890847+10*(-1+P)*P*(9071945+84*(-1+P)*P*(1310711+20*(-1+P)*P*(331975+589692*(-1+P)*P))))))))
    if(18==m):
        return -n*(-1+P)*P*(1+17*(-1+P)*P*(2027025*n**8*(-1+P)**7*P**7-8108100*n**7*(-1+P)**6*P**6*(17+74*(-1+P)*P)+630630*n**6*(-1+P)**5*P**5*(1669+4*(-1+P)*P*(4454+11347*(-1+P)*P))-10010*n**5*(-1+P)**4*P**4*(195745+2*(-1+P)*P*(2048847+16*(-1+P)*P*(808350+1586141*(-1+P)*P)))+15015*n**4*(-1+P)**3*P**3*(68713+4*(-1+P)*P*(692483+2*(-1+P)*P*(4333535+6*(-1+P)*P*(3542868+5959501*(-1+P)*P))))-26*n**3*(-1+P)**2*P**2*(5469635+2*(-1+P)*P*(235505519+2*(-1+P)*P*(2871608777+18*(-1+P)*P*(1566790981+4*(-1+P)*P*(1660236711+2531171618*(-1+P)*P)))))+3*n**2*(-1+P)*P*(1240383+4*(-1+P)*P*(75732128+3*(-1+P)*P*(1333746277+4*(-1+P)*P*(6454592891+48*(-1+P)*P*(1151761142+45*(-1+P)*P*(99061346+141547813*(-1+P)*P))))))+30*(257+4*(-1+P)*P*(63047+3276*(-1+P)*P*(2455+132*(-1+P)*P*(713+20*(-1+P)*P*(569+252*(-1+P)*P*(17+20*(-1+P)*P*(3+4*(-1+P)*P)))))))+n*(-7709-2*(-1+P)*P*(5643391+36*(-1+P)*P*(24051175+2*(-1+P)*P*(512494351+20*(-1+P)*P*(437838307+252*(-1+P)*P*(13716121+20*(-1+P)*P*(2506191+3431678*(-1+P)*P)))))))))
   
        
