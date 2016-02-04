#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include <stdlib.h> //srand, rand

void printdevice();
void initialize_neutrons(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float, int banksize,int num_src, int seed);
void initialize_neutrons_active(MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
void initialize_neutrons_active_not_src(unsigned gridx, unsigned blockx,MemStruct DeviceMem, int num_seg, int seed);
void initialize_neutrons_fix(unsigned gridx, unsigned blockx,MemStruct DeviceMem,float width,int ubat);
#if defined(__SCATTERPLOT)
void copyinitial(MemStruct DeviceMem, MemStruct HostMem, unsigned gridsize);
#endif


unsigned start_neutrons_f2(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin, MemStruct HostMem);
void start_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat,unsigned num_src,unsigned banksize, unsigned tnum_bin);
void start_neutrons_active(unsigned ibat, unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned num_seg, unsigned banksize, unsigned tnum_bin, MemStruct HostMem);

unsigned count_neutrons(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem, unsigned num_src);
unsigned count_lives(unsigned gridx, unsigned blockx, MemStruct DeviceMem, MemStruct HostMem);

void save_results(unsigned ibat, unsigned gridx, unsigned blockx, unsigned num_bin, MemStruct DeviceMem, MemStruct HostMem);
void print_results(unsigned meshes, unsigned nbat, double *tally);

unsigned setbank_converge(MemStruct DeviceMem, MemStruct HostMem, int gridsize);
void setbank_prepbank(MemStruct DeviceMem, MemStruct HostMem, int gridsize, unsigned ibat);

void setbank_active_in(unsigned idelta, MemStruct DeviceMem, MemStruct HostMem, int gridsize, int shift, int num_src);
void setbank_active_out(MemStruct DeviceMem, MemStruct HostMem, unsigned num_src, unsigned idelta);
void set_cursor_safe(MemStruct HostMem, unsigned ibat);
void bank_print(MemStruct HostMem);


void check(unsigned gridx, unsigned blockx, MemStruct DeviceMem, unsigned ubat);

int flushbank(MemStruct DeviceMem, MemStruct HostMem, unsigned lastpop,float a,unsigned gridsize,int ibat,int nbat);
int count_pop(int *live, int gridsize);
void resetcount(MemStruct DeviceMem);
#endif
