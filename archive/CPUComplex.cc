/*
================================================================================
Description: complex number class
History****:*Date*******************Author**************************************
*****Thu Apr 12 20:02:30 2012 +0800 Jilang Miao<wsmjl@mail.ustc.edu.cn> 
  A whopper with all thought out: 
  complex object with real, imag, norm and arg;
  constructor with (real,imag) or (norm,arg);
  division with 0-divisor treatment.   
*****Fri Mar 28 09:10:06 2014 -0400 Jilang Miao<jlmiao@mit.edu>
  Simplified for specific use: 
  cuda host and device; 
  no complex division needed;
  transfer norm and arg use to Norm(), Arg() functions.
*****Sat Apr 19 22:07:11 2014 -0400 Jilang Miao 缪佶朗 <jlmiao@mit.edu>
  Tailor to be compatible with general functions:
  real(complex),imag(complex),exp(complex);
  overload negative and minus operator simultaneously
================================================================================
*/
#include "CPUComplex.h"
CPUComplex::CPUComplex(double real, double imag){
  m_real=real;
  m_imag=imag;
}


CPUComplex::CPUComplex(NormArg NA){
  m_real= NA.norm * cos(NA.arg);
  m_imag=NA.norm * sin(NA.arg);
}

double CPUComplex::Re(){
  return m_real;
}


double real(CPUComplex c){
  return c.Re();
}

double CPUComplex::Im(){
  return m_imag;
}



double imag(CPUComplex c){
  return c.Im();
}

double CPUComplex::Norm(){
  return sqrt( m_real*m_real + m_imag*m_imag);
}

double Norm(CPUComplex c){
  double x = c.Re();
  double y = c.Im();
  return sqrt( x*x + y*y);
}

double CPUComplex::Arg(){
  if ( 0== m_real){
    if( m_imag>0) return PI*0.5;
    else if (m_imag<0) return PI*1.5;
    else return 0.0;
  }
  else{
    if(0==m_imag){
      return m_real>0 ? 0.0 : PI;
    }
    else{
      if ( m_real <0) return atan(m_imag/m_real )+PI;
      else return m_imag>0 ? atan(m_imag/m_real) : atan(m_imag/m_real) +2.0*PI;
    }
  }
}


double Arg(CPUComplex c){
  double x = c.Re();
  double y = c.Im();
  if ( 0 == x){
    if( y>0) return PI*0.5;
    else if (y<0) return PI*1.5;
    else return 0.0;
  }
  else{
    if(0==y){
      return x>0 ? 0.0 : PI;
    }
    else{
      if ( x <0) return atan(y/x )+PI;
      else return y>0 ? atan(y/x) : atan(y/x) +2.0*PI;
    }
  }
}

CPUComplex CPUComplex::Conjugate(){
  return CPUComplex(m_real, 0.0 - m_imag);
}

CPUComplex Conjugate(CPUComplex c){
  return CPUComplex(c.Re(), 0.0 - c.Im());
}


CPUComplex exp(CPUComplex c){
  double y = c.Im();
  return exp(c.Re())*(cos(y)*sin(y));
}

void CPUComplex::output(){
  if(0==m_real){
    if(0==m_imag) cout<<m_real<<endl;
    else cout<<m_imag<<"i"<<endl;
  }
  else{
    if(0==m_imag) cout<<m_real<<endl;
    else{
      if(m_imag>0)   cout<<m_real<<"+"<<m_imag<<"i"<<endl;
      else cout<<m_real<<"-"<<0-m_imag<<"i"<<endl;
    }
  }
}

void output(CPUComplex c){
  if(0==c.Re()){
    if(0==c.Im()) cout<<c.Re()<<endl;
    else cout<<c.Im()<<"i"<<endl;
  }
  else{
    if(0==c.Im()) cout<<c.Re()<<endl;
    else{
      if(c.Im()>0)   cout<<c.Re()<<"+"<<c.Im()<<"i"<<endl;
      else cout<<c.Re()<<"-"<<0-c.Im()<<"i"<<endl;
    }
  }
}

void CPUComplex::display(){
  double norm = Norm();
  double arg  = Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

void display(CPUComplex c){
  double norm = c.Norm();
  double arg  = c.Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

CPUComplex CPUComplex::operator +(CPUComplex c){
   return CPUComplex( m_real+c.Re(), m_imag+c.Im());
}

CPUComplex CPUComplex::operator +(double d){
  return CPUComplex(m_real+d, m_imag);
}

CPUComplex CPUComplex::operator - (double d) {
  return CPUComplex(m_real-d, m_imag);
}

 CPUComplex CPUComplex ::operator * (CPUComplex c) {
   double real=m_real*c.Re()-m_imag*c.Im();
   double imag=c.Re()*m_imag+m_real*c.Im();
  return CPUComplex(real, imag);
}

 CPUComplex CPUComplex :: operator * (double d){
  return CPUComplex( m_real*d, m_imag*d);
}


 CPUComplex CPUComplex::operator / (double d) {
  return CPUComplex( m_real/ d, m_imag/d);
}

 CPUComplex CPUComplex ::operator /(CPUComplex c) {
  double r=c.Norm() * c.Norm();
  CPUComplex c1(m_real, m_imag);
  CPUComplex c2(c.Re(), 0.0 - c.Im());
  CPUComplex c3=c1*c2;
  return c3/r;

}


/*
CPUComplex & CPUComplex ::operator =(const CPUComplex &c_right){
	if( this==&c_right) return *this;
	m_real=c_right.m_real;
	m_imag=c_right.m_imag;
	m_NA=c_right.m_NA;

	return *this;
}*/


 CPUComplex::~CPUComplex(){

}



CPUComplex operator + (double d, CPUComplex c){
  return CPUComplex( d+c.Re(), c.Im());
}

CPUComplex operator - (double d, CPUComplex c) {
  return CPUComplex (d-c.Re(), -c.Im());
}

CPUComplex operator - (CPUComplex c1, CPUComplex c2) {
  return CPUComplex (c1.Re()-c2.Re(), c1.Im()-c2.Im());
}

CPUComplex operator - (CPUComplex c){
  return CPUComplex(-c.Re(), -c.Im());
}


CPUComplex operator * (double d, CPUComplex c){
  return CPUComplex (d*c.Re(), d*c.Im());
}

 CPUComplex operator / (double d, CPUComplex c){
  return d*c.Conjugate() / (c.Norm() *c.Norm() );
}

