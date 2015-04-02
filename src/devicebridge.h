#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float);

void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned devstep);

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);

void print_results(unsigned gridx, unsigned blockx, unsigned num_src, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem, float timems);
#endif
