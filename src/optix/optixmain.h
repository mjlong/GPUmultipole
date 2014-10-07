#ifndef __OPTIXMAIN_H__
#define __OPTIXMAIN_H__
#include <optix_cuda_interop.h>

void createContext( int width, float R1, float Hh,unsigned int num, RTcontext* context, RTbuffer* buffer, RTbuffer*, RTbuffer*, float*);
void createMaterial( RTcontext context, RTmaterial* material);
void createGeometryBox( RTcontext context, RTgeometry* box, float* );
void createGeometryCylinder( RTcontext context, RTgeometry* cylinder, float,float,float,float,float,float,float);
void createGeometrySphere( RTcontext context, RTgeometry* sphere, float*);
void createInstances( RTcontext context, RTmaterial material, float*, int n, int m  );
void printUsageAndExit( const char* argv0 );
void tracemain(int, int, int, float*);

#endif
