#ifndef __QUICKW_H__
#define __QUICKW_H__

#if defined(__CFLOAT)
#define CMPTYPE float
#else
#define CMPTYPE double 
#endif

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

#define LENGTH 62
#define WIDTH  0.1

//__device__ void initialize_w_tabulated(CComplex*);
__device__ void fill_w_tabulated(CComplex<CMPTYPE>*, int);
#if defined(__QUICKWG)
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE>, CComplex<CMPTYPE>*);
#else //__QUICKWT or __QUICKWC
__device__ CComplex<CMPTYPE> w_function(CComplex<CMPTYPE>);
#endif




#endif
