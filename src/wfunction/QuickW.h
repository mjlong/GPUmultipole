#ifndef __QUICKW_H__
#define __QUICKW_H__

#if defined(__CFLOAT)
#define CMPTYPE float
#define CMPTYPE2 float2
#else
#define CMPTYPE double 
#define CMPTYPE2 double2
#endif

#include <stdio.h>
#include "global.h"
#include "CComplex.h"
#include "Faddeeva.h"

/*=============================================================================== 
  FOR ALL FUNCTIONS BELOW:             
  Copyright (c) 2011-2014 Massachusetts Institute of Technology            
  Permission is hereby granted, free of charge, to any person obtaining     
  a copy of this software and associated documentation files (the                 
  "Software"), to deal in the Software without restriction, including                  
  without limitation the rights to use, copy, modify, merge, publish,             
  distribute, sublicense, and/or sell copies of the Software, and to              
  permit persons to whom the Software is furnished to do so, subject to                           
  the following conditions:                             
                                                                                                   
  The above copyright notice and this permission notice shall be                 
  included in all copies or substantial portions of the Software.               

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,             
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF              
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                           
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE          
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION          
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION           
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                
===============================================================================*/

#define RADIUS 6.0
#define WIDTH  0.1
#define LENGTH 62

//__device__ void initialize_w_tabulated(CComplex*);
__global__ void fill_w_tabulated(CComplex<CMPTYPE>*);

#if defined(__QUICKWG)
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE>, CComplex<CMPTYPE>*);
#elif defined(__QUICKWT)
__host__ void bindwtable(CComplex<CMPTYPE>* wtable);
__host__ void unbindwtable();
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE>);
#else
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE> );
#endif




#endif
