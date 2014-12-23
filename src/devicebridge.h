#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
#include "multipole_data.h"
#include "multipole.h"
#include "material_data.h"
#include "material.h"

void printdevice();
void print_results(unsigned num_src, unsigned num_bin, MemStruct HostMem, float timems);
#endif
