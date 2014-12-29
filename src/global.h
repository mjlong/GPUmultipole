#if defined(__ALLCPU)
#include "CPUComplex.h"
#define ccomplex CPUComplex<CMPTYPE>
#endif
#if defined(__XS_GPU)||defined(__W__GPU)||defined(__PFOURIERW)
//#include "CComplex.h"
#define ccomplex CComplex<CMPTYPE>
#endif


#define KB (8.617342E-5)
#define C1 0.002196771 
#define C2 1301997.933
#define ONETRD  1.0/3.0
#define HALFI ccomplex(0.0,0.5)
#define ONEI  ccomplex(0.0,1.0)


#if defined(__CFLOAT)
#define CMPTYPE float
#else
#define CMPTYPE double
#endif
