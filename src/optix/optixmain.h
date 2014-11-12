#ifndef __OPTIXMAIN_H__
#define __OPTIXMAIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <optix.h>
#include "sutil.h"
#include <vector_types.h>
#include <math.h>
#include "commonStructs.h"
#include <optix_cuda_interop.h>

#include "neutron.h"
#include "multipole.h"

void createMaterial( RTcontext context, RTmaterial* material);
void createGeometryBox( RTcontext context, RTgeometry* box, float* );
void createGeometryCylinder( RTcontext context, RTgeometry* cylinder, float,float,float,float,float,float,float);
void createGeometrySphere( RTcontext context, RTgeometry* sphere, float*);
void createInstances( RTcontext context, RTmaterial material, float*, int n, int m  );
void tracemain(int, int, int, float*, NeutronInfoStruct);
#if defined(__QUICKW)
void initialize_context(RTcontext, int, int, int, float*, NeutronInfoStruct,multipole, CComplex<double>*);
void createContext( int width, float R1, float Hh, unsigned num_geo, RTcontext context, NeutronInfoStruct nInfo,multipole, CComplex<double>*);
#else
void initialize_context(RTcontext, int, int, int, float*, NeutronInfoStruct,multipole);
void createContext( int width, float R1, float Hh, unsigned num_geo, RTcontext context, NeutronInfoStruct nInfo,multipole);
#endif


#endif
