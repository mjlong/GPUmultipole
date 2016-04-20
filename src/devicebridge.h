#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include <stdlib.h> //srand, rand

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float, int banksize,int num_src, int seed);
void initialize_neutrons_fix(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int ubat);
#if defined(__SCATTERPLOT)
void copyinitial(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize);
#endif


#if defined(__FTALLY2)
unsigned start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin, MemStruct HostMem);
#else
void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin);
#endif


unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);

void save_results(unsigned ibat, unsigned gridx, unsigned blockx, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem);
void print_results(unsigned meshes, unsigned nbat, double *tally);

#if defined(__MTALLY)||(__FTALLY)||(__FTALLY2)
unsigned setbank_converge(MemStruct DeviceMem, MemStruct HostMem, int gridsize);
void copysrcforwrite(MemStruct HostMem, int num_src, float* x2, float* y2, float* z2);
#if defined(FTALLY2)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize,
		 int banksize, unsigned jstart, int shift);
void setbank2(MemStruct DeviceMem, MemStruct HostMem, int banksize, unsigned jstart);
#else// MTALLY or FTALLY
#if defined(__MTALLY)
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize,
		 int oldbanksize, int tnum_bins);
#else
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, int tnum_bins);
#endif
#endif
#else//CTALLY
unsigned setbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize);
#endif
void check(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat);

int flushbank(MemStruct DeviceMem, MemStruct HostMem, unsigned lastpop,float a,unsigned gridsize,int ibat,int nbat);
int count_pop(int *live, int gridsize);
void resetcount(MemStruct DeviceMem);
#endif
