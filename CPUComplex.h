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
#ifndef __CPUCOMPLEX_H__
#define __CPUCOMPLEX_H__
#include <iostream>
#include <math.h>
using namespace std;

#define PI 3.14159265358979323846264338327950


struct NormArg{
  double norm;
  double arg;
};

class CPUComplex{
private:
  double m_real;
  double m_imag;
  //  NormArg m_NA;

public:
   CPUComplex(double m_real=0.0, double m_imag=0.0f);
   CPUComplex(NormArg m_NA );
   double Re();
   double Im();
   double Norm();
   double Arg();
   CPUComplex Conjugate();
   //display() in form of norm*Exp(phi)
   void display();
   //output()  in form of real+(-)imag(*)i
   void output();
   ~CPUComplex();
   CPUComplex operator + (CPUComplex);
   CPUComplex operator + (double);
   CPUComplex operator -  (double);
   CPUComplex operator *  (CPUComplex);
   CPUComplex operator *  (double);
  //division operator has been simplified by ignoring 0 divisor
  //in case of being devided by double, it is suggested to multiply by 1/double
   CPUComplex operator /   (CPUComplex);
   CPUComplex operator /   (double);
  //	CPUComplex & operator =  (const CPUComplex & c);

};

void output(CPUComplex);
void display(CPUComplex);
CPUComplex Conjugate(CPUComplex);
double Norm(CPUComplex);
double Arg(CPUComplex);
double real(CPUComplex);
double imag(CPUComplex);
CPUComplex exp(CPUComplex);

CPUComplex operator + (double d, CPUComplex c);
CPUComplex operator - (double d, CPUComplex c);
CPUComplex operator -  (CPUComplex);
//It should be noted that the '-' here is negative sign
CPUComplex operator - (CPUComplex c1, CPUComplex c2);
CPUComplex operator * (double d, CPUComplex c);
CPUComplex operator / (double d, CPUComplex c);


#endif
