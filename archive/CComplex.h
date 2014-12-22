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
#ifndef __CCOMPLEX_H__
#define __CCOMPLEX_H__
#include <iostream>
#include <math.h>
using namespace std;

#define PI 3.141592653589793238462643383279


struct NormArg{
  double norm;
  double arg;
};

class CComplex{
private:
  double m_real;
  double m_imag;
  //  NormArg m_NA;

public:
  __host__ __device__ CComplex(double m_real=0.0, double m_imag=0.0f);
  __host__ __device__ CComplex(NormArg m_NA );
  __host__ __device__ double Re();
  __host__ __device__ double Re(CComplex);
  __host__ __device__ double Im();
  __host__ __device__ double Im(CComplex);
  __host__ __device__ double Norm();
  __host__ __device__ double Norm(CComplex);
  __host__ __device__ double Arg();
  __host__ __device__ double Arg(CComplex);
  __host__ __device__ CComplex Conjugate();
  __host__ __device__ CComplex Conjugate(CComplex);
  //display() in form of norm*Exp(phi)
  void display();
  void display(CComplex);
  //output()  in form of real+(-)imag(*)i
  void output();
  void output(CComplex);
  __host__ __device__ ~CComplex();
  __host__ __device__ CComplex operator + (CComplex);
  __host__ __device__ CComplex operator + (double);
  __host__ __device__ CComplex operator -  (CComplex);
  __host__ __device__ CComplex operator -  (double);
  __host__ __device__ CComplex operator *  (CComplex);
  __host__ __device__ CComplex operator *  (double);
  //division operator has been simplified by ignoring 0 divisor
  //in case of being devided by double, it is suggested to multiply by 1/double
  __host__ __device__ CComplex operator /   (CComplex);
  __host__ __device__ CComplex operator /   (double);
  //	CComplex & operator =  (const CComplex & c);

};

__host__ __device__ CComplex operator + (double d, CComplex c);
__host__ __device__ CComplex operator - (double d, CComplex c);
__host__ __device__ CComplex operator * (double d, CComplex c);
__host__ __device__ CComplex operator / (double d, CComplex c);


#endif
