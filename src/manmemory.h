#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "gpuerrchk.h"
#include "neutron.h"
#include "CComplex.h"

//initialize device
void initialize_device();

//Simulation memory allocate and deallocate
void initialize_memory(MemStruct *HostMem, unsigned numbins);
void release_memory(MemStruct DeviceMem);

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

#if defined(__XS_GPU)
void allocate_buffer(unsigned maxiso, unsigned** iS_d, 
                     CMPTYPE** sigTs_h, CMPTYPE** sigAs_h, CMPTYPE** sigFs_h,
                     CMPTYPE** sigTs_d, CMPTYPE** sigAs_d, CMPTYPE** sigFs_d);
void release_buffer(unsigned* iS_d, 
                    CMPTYPE* sigTs_h, CMPTYPE* sigAs_h, CMPTYPE* sigFs_h, 
                    CMPTYPE* sigTs_d, CMPTYPE* sigAs_d, CMPTYPE* sigFs_d);
#endif

#endif
