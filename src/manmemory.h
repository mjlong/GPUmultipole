#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "gpuerrchk.h"
#include "neutron.h"

//initialize device
void initialize_device();

void initialize_memory_bank(MemStruct *HostMem, unsigned banksize);
void release_memory_bank(MemStruct HostMem);
//Simulation memory allocate and deallocate
void initialize_memory_data(MemStruct *DeviceMem, MemStruct *HostMem);
void copydata(MemStruct DeviceMem, MemStruct HostMem);
void release_memory_data(MemStruct DeviceMem, MemStruct HostMem);
void allocate_memory_converge(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned num_seg);
void   allocate_memory_active(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx,unsigned num_seg);
void release_memory_converge(MemStruct DeviceMem, MemStruct HostMem);
void   release_memory_active(MemStruct DeviceMem, MemStruct HostMem);

void initialize_memory(MemStruct *DeviceMem, MemStruct *HostMem, unsigned numbins, unsigned gridx, unsigned blockx, unsigned numbatches,unsigned ubat);
void release_memory(MemStruct DeviceMem, MemStruct HostMem);
void copymeans(int *h_cnt, int *batcnt, unsigned meshes, unsigned offset);
void resettally(CMPTYPE *cnt, unsigned totbins);
void resettally(int *cnt, unsigned totbins);
#endif
