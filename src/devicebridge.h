#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include "multipole_data.h"
#include "multipole.h"
#include "material_data.h"
#include "material.h"

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem);

void start_neutrons(unsigned gridx, unsigned blockx, unsigned numIsos, multipole mp_para, CMPTYPE* devicearray, MemStruct DeviceMem, unsigned num_src, unsigned devstep);

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);

void remain_neutrons(unsigned gridx, unsigned blockx, unsigned numIsos, multipole mp_para, CMPTYPE* devicearray, MemStruct DeviceMem);

void print_results(unsigned gridx, unsigned blockx, unsigned num_src, unsigned devstep, MemStruct DeviceMem, MemStruct HostMem, CMPTYPE* hostarray, CMPTYPE* devicearray, unsigned* blockcnt,unsigned* cnt, float timems);


#endif
