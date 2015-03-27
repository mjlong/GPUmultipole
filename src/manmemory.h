#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "gpuerrchk.h"
#include "neutron.h"
#include "CComplex.h"
#include "CPUComplex.h"

//initialize device
void initialize_device();

void copyE(MemStruct HostMem, MemStruct DeviceMem, unsigned gridsize);
void copyZ(CPUComplex<CMPTYPE>* pz_h, CComplex<CMPTYPE>* pz_d, unsigned numz);
void copyW(CComplex<CMPTYPE>* pw_d, CPUComplex<CMPTYPE>* pw_h, unsigned numz);
void allocateZW(CComplex<CMPTYPE> **pz, CComplex<CMPTYPE> **pw, unsigned numz);
void releaseZW(CComplex<CMPTYPE>* pz, CComplex<CMPTYPE>* pw);

//Simulation memory allocate and deallocate
void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx );
void release_memory(MemStruct DeviceMem, MemStruct HostMem);

//Faddeeva function table management
#if defined(__FOURIERW)
void fill_wtables(CMPTYPE** da, CMPTYPE** db);
void release_wtables(CMPTYPE* da, CMPTYPE* db);
#endif

#if defined(__INTERPEXP)
void fill_wtables(CComplex<CMPTYPE>** exptable);
void release_wtables(CComplex<CMPTYPE>* exptable);
#endif

#if defined(__QUICKW)
void fill_wtables(CComplex<CMPTYPE>** wtable);
void release_wtables(CComplex<CMPTYPE>* wtable);
#endif



#endif
