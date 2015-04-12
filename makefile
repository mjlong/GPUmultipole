DIR_SRC = ./src
DIR_SRC_PTX = ./src/ptx
DIR_SRC_OPT = ./src/optix
DIR_OBJ = ./obj
DIR_PTX = ./obj/ptx
DIR_HDF5  = /opt/hdf5/1.8.14-gnu#/home/jlmiao/opt/hdf5
DIR_CUDA6 = /home/jmiao/cuda-6.0
DIR_OPTIX = /home/jlmiao/Documents/NVIDIA-OptiX-SDK-3.6.0-linux64
#Include flags
INC_SRC   = -I${DIR_SRC} 
INC_HDF5  = -I${DIR_HDF5}/include
INC_CUDA6 = -I${DIR_CUDA6}/include
INC_OPTIX = -I${DIR_OPTIX}/include -I${DIR_SRC_OPT}
NCINCFLAGS  = $(INC_SRC) $(INC_CUDA6) $(INC_OPTIX)
CCINCFLAGS  = $(INC_SRC) $(INC_HDF5) $(INC_OPTIX)
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
ifeq ($(ver),debug)
NCFLAGS=-g -G -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
CCFLAGS=-c -g                 $(CCINCFLAGS) 
DIR_BIN = ./bin/debug
else
NCFLAGS=      -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
CCFLAGS=-c                    $(CCINCFLAGS) 
DIR_BIN = ./bin/release
endif
LINKLAG=   -dlink -arch=sm_20  
LDFLAGS=-L${DIR_HDF5}/lib/ -L${DIR_CUDA6}/lib64 -L${DIR_OPTIX}/lib64/ -loptix -lcudart -lhdf5 -lstdc++
GSOURCES=$(wildcard ${DIR_SRC}/*.cu)
PSOURCES=$(wildcard ${DIR_SRC_PTX}/*.cu)
ifeq ($(version),many)
RTMETHOD = -D __MANY__
PTXFIX =_many.ptx
EXEFIX = _many
else
RTMETHOD=
PTXFIX =_one.ptx
EXEFIX = _one
endif
EXENAME=$(DIR_BIN)/gpu_box
EXECUTABLE=$(EXENAME)

CSOURCES=$(wildcard ${DIR_SRC}/*.cc)
MSOURCES=$(wildcard ${DIR_SRC}/*.cxx)
CNVCCSRC=$(wildcard ${DIR_SRC_OPT}/*.cxx)
COBJECTS=$(patsubst %.cc, ${DIR_OBJ}/%.obj, $(notdir ${CSOURCES}))
MOBJECTS=$(patsubst %.cxx, ${DIR_OBJ}/%.ob, $(notdir ${MSOURCES}))
GOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${GSOURCES}))
PTXJECTS=$(patsubst ${DIR_SRC_PTX}/%.cu, ${DIR_PTX}/%$(PTXFIX),  ${PSOURCES})
CNVCCOBJ=$(patsubst %.cxx, ${DIR_OBJ}/%.ob, $(notdir ${CNVCCSRC}))
LINKJECT=${DIR_OBJ}/dlink.o      
all: $(PTXJECTS) $(EXECUTABLE)
${DIR_PTX}/%$(PTXFIX) : ${DIR_SRC_PTX}/%.cu
	$(NVCC) $(RTMETHOD) $(NCPLAGS) -use_fast_math $^ -o $@
$(EXECUTABLE): $(MOBJECTS) $(COBJECTS) $(CNVCCOBJ) $(GOBJECTS) $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@
${DIR_OBJ}/%.obj : ${DIR_SRC}/%.cc
	$(CC)             $(CMPTYPE) $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC}/%.cxx
	$(NVCC) $(RTMETHOD) $(CMPTYPE) $(NCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	$(NVCC) $(RTMETHOD) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/wfunction/%.cu
	$(NVCC) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC_OPT}/%.cxx
	$(NVCC) $(RTMETHOD) $(NCFLAGS) $^ -o $@
$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.ob  -exec rm -rf {} \;
clean :  
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.ob  -exec rm -rf {} \;
	find ${DIR_PTX} -name *.ptx -exec rm -rf {} \;
