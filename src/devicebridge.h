#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include "multipole_data.h"
#include "multipole.h"
#include "material_data.h"
#include "material.h"

#include "Faddeeva.hh"

void z2w(CPUComplex<CMPTYPE> *pz, CPUComplex<CMPTYPE>* pw, unsigned numz);
void z2w(CComplex<CMPTYPE> *pz,   CComplex<CMPTYPE>* pw,   unsigned numz);
void call_idle_kernel();

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem);

void start_neutrons(unsigned gridx, unsigned blockx, material mat, multipole mp_para, MemStruct DeviceMem, unsigned num_src,unsigned active);

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);
void sort_prepare(unsigned gridx, unsigned blockx,MemStruct DeviceMem, material mat);

void transport_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem, material mat,unsigned renew);

void print_results(unsigned gridx, unsigned blockx, unsigned num_src, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem, float timems);
#endif
