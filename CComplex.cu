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
================================================================================
*/
#include "CComplex.h"
 __host__ __device__ CComplex::CComplex(double real, double imag){
  m_real=real;
  m_imag=imag;
}


 __host__ __device__ CComplex::CComplex(NormArg NA){
  m_real= NA.norm * cos(NA.arg);
  m_imag=NA.norm * sin(NA.arg);
}

__host__ __device__ double CComplex::Re(){
  return m_real;
}

__host__ __device__ double CComplex::Re(CComplex c){
  return c.m_real;
}

__host__ __device__ double CComplex::Im(){
  return m_imag;
}

__host__ __device__ double CComplex::Im(CComplex c){
  return c.m_imag;
}

__host__ __device__ double CComplex::Norm(){
  return sqrt( m_real*m_real + m_imag*m_imag);
}

__host__ __device__ double CComplex::Norm(CComplex c){
  return sqrt( c.m_real*c.m_real + c.m_imag*c.m_imag);
}

__host__ __device__ double CComplex::Arg(){
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


__host__ __device__ double CComplex::Arg(CComplex c){
  if ( 0 == c.m_real){
    if( c.m_imag>0) return PI*0.5;
    else if (c.m_imag<0) return PI*1.5;
    else return 0.0;
  }
  else{
    if(0==c.m_imag){
      return c.m_real>0 ? 0.0 : PI;
    }
    else{
      if ( c.m_real <0) return atan(c.m_imag/c.m_real )+PI;
      else return c.m_imag>0 ? atan(c.m_imag/c.m_real) : atan(c.m_imag/c.m_real) +2.0*PI;
    }
  }
}

__host__ __device__ CComplex CComplex::Conjugate(){
  return CComplex(m_real, 0.0 - m_imag);
}

__host__ __device__ CComplex CComplex::Conjugate(CComplex c){
  return CComplex(c.m_real, 0.0 - c.m_real);
}

void CComplex::output(){
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

void CComplex::output(CComplex c){
  if(0==c.m_real){
    if(0==c.m_imag) cout<<c.m_real<<endl;
    else cout<<c.m_imag<<"i"<<endl;
  }
  else{
    if(0==c.m_imag) cout<<c.m_real<<endl;
    else{
      if(c.m_imag>0)   cout<<c.m_real<<"+"<<c.m_imag<<"i"<<endl;
      else cout<<c.m_real<<"-"<<0-c.m_imag<<"i"<<endl;
    }
  }
}

void CComplex::display(){
  double norm = Norm();
  double arg  = Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

void CComplex::display(CComplex c){
  double norm = c.Norm();
  double arg  = c.Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

__host__ __device__ CComplex CComplex::operator +(CComplex c){
  return CComplex( m_real+c.m_real, m_imag+c.m_imag);
}

__host__ __device__ CComplex CComplex::operator +(double d){
  return CComplex(m_real+d, m_imag);
}

__host__ __device__ CComplex CComplex::operator - (CComplex c){
  return CComplex(m_real-c.m_real, m_imag-c.m_imag);
}

__host__ __device__ CComplex CComplex::operator - (double d) {
  return CComplex(m_real-d, m_imag);
}

__host__ __device__ CComplex CComplex ::operator * (CComplex c) {
  double real=m_real*c.m_real-m_imag*c.m_imag;
  double imag=c.m_real*m_imag+m_real*c.m_imag;
  return CComplex(real, imag);
}

__host__ __device__ CComplex CComplex :: operator * (double d){
  return CComplex( m_real*d, m_imag*d);
}


__host__ __device__ CComplex CComplex::operator / (double d) {
  return CComplex( m_real/ d, m_imag/d);
}

__host__ __device__ CComplex CComplex ::operator /(CComplex c) {
  double r=c.Norm() * c.Norm();
  CComplex c1(m_real, m_imag);
  CComplex c2(c.m_real, 0.0 - c.m_imag);
  CComplex c3=c1*c2;
  return c3/r;

}


/*
CComplex & CComplex ::operator =(const CComplex &c_right){
	if( this==&c_right) return *this;
	m_real=c_right.m_real;
	m_imag=c_right.m_imag;
	m_NA=c_right.m_NA;

	return *this;
}*/


__host__ __device__ CComplex::~CComplex(){

}



__host__ __device__ CComplex operator + (double d, CComplex c){
  return CComplex( d+c.Re(), c.Im());
}

__host__ __device__ CComplex operator - (double d, CComplex c) {
  return CComplex (d-c.Re(), c.Im());
}

__host__ __device__ CComplex operator * (double d, CComplex c){
  return CComplex (d*c.Re(), d*c.Im());
}

__host__ __device__ CComplex operator / (double d, CComplex c){
  return d*c.Conjugate() / (c.Norm() *c.Norm() );
}

