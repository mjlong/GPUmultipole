#ifndef __DEVICEBRIDGE_H__
#define __DEVICEBRIDGE_H__
void printdevice();
void uploadmultipole(struct multipoledata* data, unsigned numIsos);

void anyvalue(struct multipoledata* data, unsigned numIsos, struct matdata* pmat, unsigned totIsos, unsigned setgridx, unsigned setblockx, unsigned num_src, unsigned devstep, unsigned* cnt, unsigned* blockcnt, CMPTYPE* hostarray, CMPTYPE* devicearray, MemStruct HostMem, MemStruct DeviceMem);

#endif
