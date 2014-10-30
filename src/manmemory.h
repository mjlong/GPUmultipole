#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "gpuerrchk.h"
#include "neutron.h"

//initialize device
void initialize_device();

//Simulation memory allocate and deallocate
void initialize_memory(MemStruct* DeviceMem, MemStruct* HostMem, CMPTYPE** devicearray, CMPTYPE** hostarray, unsigned **cnt, unsigned** blockcnt, unsigned gridx, unsigned blockx );

void release_memory(MemStruct DeviceMem, MemStruct HostMem, CMPTYPE* devicearray, CMPTYPE* hostarray, unsigned *cnt, unsigned* blockcnt );

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
