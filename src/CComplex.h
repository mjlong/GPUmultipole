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
*****Thu Jun 5 22:47:23 2014 -0400  Jilang Miao 缪佶朗<jlmiao@mit.edu>
  Make it a class template
================================================================================
*/
#ifndef __CCOMPLEX_H__
#define __CCOMPLEX_H__
#include <iostream>
#include <math_functions_dbl_ptx3.h>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846264338327950
#endif

template <class DType>
class CComplex{
private:
  DType m_real;
  DType m_imag;

public:
   __host__ __device__ CComplex(DType m_real=0.0, DType m_imag=0.0f);
   __host__ __device__ DType Re();
   __host__ __device__ DType Im();
   __host__ __device__ DType Norm();
   __host__ __device__ DType Arg();
   __host__ __device__ CComplex Conjugate();
   //display() in form of norm*Exp(phi)
   __host__ __device__ void display();
   //output()  in form of real+(-)imag(*)i
   __host__ __device__ void output();
   __host__ __device__ ~CComplex();
   __host__ __device__ CComplex operator + (CComplex);
   __host__ __device__ CComplex operator + (DType);
   __host__ __device__ CComplex operator -  (DType);
   __host__ __device__ CComplex operator *  (CComplex);
   __host__ __device__ CComplex operator *  (DType);
  //division operator has been simplified by ignoring 0 divisor
  //in case of being devided by DType, it is suggested to multiply by 1/DType
   __host__ __device__ CComplex operator /   (CComplex);
   __host__ __device__ CComplex operator /   (DType);
  //	CComplex & operator =  (const CComplex & c);

};

template <class DType>
__host__ __device__ CComplex<DType>::CComplex(DType real, DType imag){
  m_real=real;
  m_imag=imag;
}


template <class DType>
__host__ __device__ DType CComplex<DType>::Re(){
  return m_real;
}

template <class DType>
__host__ __device__ DType real(CComplex<DType> c){
  return c.Re();
}

template <class DType>
__host__ __device__ DType CComplex<DType>::Im(){
  return m_imag;
}



template <class DType>
__host__ __device__ DType imag(CComplex<DType> c){
  return c.Im();
}

template <class DType>
__host__ __device__ DType CComplex<DType>::Norm(){
  return sqrt( m_real*m_real + m_imag*m_imag);
}

template <class DType>
__host__ __device__ DType Norm(CComplex<DType> c){
  DType x = c.Re();
  DType y = c.Im();
  return sqrt( x*x + y*y);
}

template <class DType>
__host__ __device__ DType CComplex<DType>::Arg(){
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


template <class DType>
__host__ __device__ DType Arg(CComplex<DType> c){
  DType x = c.Re();
  DType y = c.Im();
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

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType>::Conjugate(){
  return CComplex(m_real, 0.0 - m_imag);
}

template <class DType>
__host__ __device__ CComplex<DType> Conjugate(CComplex<DType> c){
  return CComplex<DType>(c.Re(), 0.0 - c.Im());
}


template <class DType>
__host__ __device__ CComplex<DType> exp(CComplex<DType> c){
  return CComplex<DType>(exp(c.Re())*cosf(c.Im()),exp(c.Re())*sinf(c.Im()));
}

template <class DType>
__host__ __device__ void CComplex<DType>::output(){
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

template <class DType>
__host__ __device__ void output(CComplex<DType> c){
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

template <class DType>
__host__ __device__ void CComplex<DType>::display(){
  DType norm = Norm();
  DType arg  = Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

template <class DType>
__host__ __device__ void display(CComplex<DType> c){
  DType norm = c.Norm();
  DType arg  = c.Arg();
  if(0==norm) cout<<norm<<endl;
  else cout<<norm<<"*Exp("<<arg<<")"<<endl;
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType>::operator +(CComplex<DType> c){
   return CComplex( m_real+c.Re(), m_imag+c.Im());
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType>::operator +(DType d){
  return CComplex(m_real+d, m_imag);
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType>::operator - (DType d) {
  return CComplex(m_real-d, m_imag);
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType> ::operator * (CComplex<DType> c) {
   DType real=m_real*c.Re()-m_imag*c.Im();
   DType imag=c.Re()*m_imag+m_real*c.Im();
  return CComplex(real, imag);
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType> :: operator * (DType d){
  return CComplex( m_real*d, m_imag*d);
}


template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType>::operator / (DType d) {
  return CComplex( m_real/ d, m_imag/d);
}

template <class DType>
__host__ __device__ CComplex<DType> CComplex<DType> ::operator /(CComplex<DType> c) {
  DType r=c.Norm() * c.Norm();
  CComplex<DType> c1(m_real, m_imag);
  CComplex<DType> c2(c.Re(), 0.0 - c.Im());
  CComplex<DType> c3=c1*c2;
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


template <class DType>
__host__ __device__ CComplex<DType>::~CComplex(){

}



template <class DType>
__host__ __device__ CComplex<DType> operator + (DType d, CComplex<DType> c){
  return CComplex<DType>( d+c.Re(), c.Im());
}

template <class DType>
__host__ __device__ CComplex<DType> operator - (DType d, CComplex<DType> c) {
  return CComplex<DType> (d-c.Re(), -c.Im());
}

template <class DType>
__host__ __device__ CComplex<DType> operator - (CComplex<DType> c1, CComplex<DType> c2) {
  return CComplex<DType> (c1.Re()-c2.Re(), c1.Im()-c2.Im());
}

template <class DType>
__host__ __device__ CComplex<DType> operator - (CComplex<DType> c){
  return CComplex<DType>(-c.Re(), -c.Im());
}


template <class DType>
__host__ __device__ CComplex<DType> operator * (DType d, CComplex<DType> c){
  return CComplex<DType> (d*c.Re(), d*c.Im());
}

template <class DType>
__host__ __device__ CComplex<DType> operator / (DType d, CComplex<DType> c){
  return d*c.Conjugate() / (c.Norm() *c.Norm() );
}

#endif
