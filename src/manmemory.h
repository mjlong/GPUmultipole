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
#endif
