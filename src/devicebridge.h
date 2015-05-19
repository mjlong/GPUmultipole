#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float, int banksize,int num_src);
#if defined(__SCATTERPLOT)
void copyinitial(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize);
#endif

#if defined(__TRAN)
void transient_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize,float p2);
#else
void     start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize);
#endif
unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);

void save_results(unsigned ibat, unsigned gridx, unsigned blockx, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem);
void print_results(unsigned meshes, unsigned nbat, double *tally);

unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize);

int flushbank(MemStruct DeviceMem, MemStruct HostMem, unsigned lastpop,float a,unsigned gridsize);
int count_pop(int *live, int gridsize);
void resetcount(MemStruct DeviceMem);
#endif
