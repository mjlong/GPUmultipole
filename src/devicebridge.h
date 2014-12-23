#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include "multipole_data.h"
#include "multipole.h"
#include "material_data.h"
#include "material.h"

void printdevice();
void print_results(unsigned num_src, unsigned num_bin, MemStruct HostMem, float timems);
void eval_xs(multipole mp_para, int* iS_h, int* iS_d, unsigned numIso, CMPTYPE E, CMPTYPE sqrtKT, 
                               CMPTYPE *sigTs_h, CMPTYPE *sigAs_h, CMPTYPE *sigFs_h,
                               CMPTYPE *sigTs_d, CMPTYPE *sigAs_d, CMPTYPE *sigFs_d);
#endif
