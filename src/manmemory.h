#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "global.h"
#include "gpuerrchk.h"
#include "neutron.h"

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
#if defined(__ALLCPU)
#include "CPUComplex.h"
void fill_wtables(CPUComplex<CMPTYPE>** wtable);
void release_wtables(CPUComplex<CMPTYPE>* wtable);
#else
#include "CComplex.h"
void fill_wtables(CComplex<CMPTYPE>** wtable);
void release_wtables(CComplex<CMPTYPE>* wtable);
#endif
#endif

#if defined(__XS_GPU)
void allocate_buffer(unsigned maxiso, unsigned** iS_d, 
                     CMPTYPE** sigTs_h, CMPTYPE** sigAs_h, CMPTYPE** sigFs_h,
                     CMPTYPE** sigTs_d, CMPTYPE** sigAs_d, CMPTYPE** sigFs_d);
void release_buffer(unsigned* iS_d, 
                    CMPTYPE* sigTs_h, CMPTYPE* sigAs_h, CMPTYPE* sigFs_h, 
                    CMPTYPE* sigTs_d, CMPTYPE* sigAs_d, CMPTYPE* sigFs_d);
#endif

#if defined(__W__GPU)||defined(__PFOURIERW)
#include "CComplex.h"
#include "CPUComplex.h"
void allocate_zwarray(CPUComplex<CMPTYPE>** z_h, CComplex<CMPTYPE>** z_d, 
                      CPUComplex<CMPTYPE>** w_h, CComplex<CMPTYPE>** w_d, 
                      unsigned numiso, struct multipoledata *iso);
void release_zwarray(CPUComplex<CMPTYPE>* z_h, CComplex<CMPTYPE>* z_d, 
                     CPUComplex<CMPTYPE>* w_h, CComplex<CMPTYPE>* w_d);
#endif
#endif
