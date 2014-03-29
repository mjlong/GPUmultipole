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
   CComplex(double m_real=0.0, double m_imag=0.0f);
   CComplex(NormArg m_NA );
   double Re();
   double Re(CComplex);
   double Im();
   double Im(CComplex);
   double Norm();
   double Norm(CComplex);
   double Arg();
   double Arg(CComplex);
   CComplex Conjugate();
   CComplex Conjugate(CComplex);
  //display() in form of norm*Exp(phi)
  void display();
  void display(CComplex);
  //output()  in form of real+(-)imag(*)i
  void output();
  void output(CComplex);
   ~CComplex();
   CComplex operator + (CComplex);
   CComplex operator + (double);
   CComplex operator -  (CComplex);
   CComplex operator -  (double);
   CComplex operator *  (CComplex);
   CComplex operator *  (double);
  //division operator has been simplified by ignoring 0 divisor
  //in case of being devided by double, it is suggested to multiply by 1/double
   CComplex operator /   (CComplex);
   CComplex operator /   (double);
  //	CComplex & operator =  (const CComplex & c);

};

 CComplex operator + (double d, CComplex c);
 CComplex operator - (double d, CComplex c);
 CComplex operator * (double d, CComplex c);
 CComplex operator / (double d, CComplex c);


#endif
