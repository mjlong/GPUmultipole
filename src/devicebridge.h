#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include <stdlib.h> //srand, rand

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float, int banksize,int num_src, int seed);
void initialize_neutrons_fix(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int ubat);
#if defined(__SCATTERPLOT)
void copyinitial(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize);
#endif

void     start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_src,unsigned active,unsigned banksize, unsigned tnum_bin);
unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);

void save_results(unsigned ibat, unsigned gridx, unsigned blockx, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem);
void print_results(unsigned meshes, unsigned nbat, double *tally);

#if defined(__MTALLY)||(__FTALLY)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, int tnum_bins);
#else
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize);
#endif
void check(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat);

int flushbank(MemStruct DeviceMem, MemStruct HostMem, unsigned lastpop,float a,unsigned gridsize,int ibat,int nbat);
int count_pop(int *live, int gridsize);
void resetcount(MemStruct DeviceMem);
#endif
