#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include "multipole_data.h"
#include "material_data.h"
#include "material.h"

void print_results(unsigned num_src, unsigned num_bin, MemStruct HostMem, float timems);

#if defined(__XS_GPU)
#include "multipole.h"
void printdevice();
void eval_xs(multipole mp_para, unsigned int* iS_h, unsigned int* iS_d, unsigned numIso, CMPTYPE E, CMPTYPE sqrtKT, 
                               CMPTYPE *sigTs_h, CMPTYPE *sigAs_h, CMPTYPE *sigFs_h,
                               CMPTYPE *sigTs_d, CMPTYPE *sigAs_d, CMPTYPE *sigFs_d);
#endif

#if defined(__W__GPU)
void eval_w(CPUComplex<CMPTYPE>* z_h, void** z_d, CPUComplex<CMPTYPE>* w_h, void** w_d,unsigned window);
#endif


#endif
