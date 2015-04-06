#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "gpuerrchk.h"
#include "neutron.h"

//initialize device
void initialize_device();

//Simulation memory allocate and deallocate
void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx, unsigned numbatches,unsigned ubat);
void release_memory(MemStruct DeviceMem, MemStruct HostMem);
void copymeans(unsigned *h_cnt, unsigned *acccnt, unsigned meshes, unsigned offset);
void copydata(MemStruct DeviceMem, MemStruct HostMem);
void resettally(unsigned *cnt, unsigned totbins);
#endif
