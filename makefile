ifeq ($(dim),3d)
CIDEN = -D __3D
else
CIDEN = -D __1D
endif

ifeq ($(tran),true)
CIDEN += -D __TRAN
endif

ifeq ($(ctally),true)
CIDEN += -D __TALLY
CIDEN += -D __CTALLY
endif

ifeq ($(mtally),true)
CIDEN += -D __TALLY
CIDEN += -D __MTALLY
endif

ifeq ($(ftally),true)
CIDEN += -D __TALLY
CIDEN += -D __FTALLY
endif

ifeq ($(scatterplt),true)
CIDEN += -D __SCATTERPLOT
endif

DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_HDF5  = /opt/hdf5/1.8.14-gnu#/home/jlmiao/opt/hdf5
DIR_CUDA6 = /home/jmiao/cuda-6.0
#Include flags
INC_SRC   = -I${DIR_SRC} 
INC_HDF5  = -I${DIR_HDF5}/include
INC_CUDA6 = -I${DIR_CUDA6}/include
NCINCFLAGS  = $(INC_SRC) $(INC_CUDA6) 
CCINCFLAGS  = $(INC_SRC) $(INC_HDF5) 
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
LDFLAGS=-L${DIR_HDF5}/lib/ -L${DIR_CUDA6}/lib64  -lcudart -lhdf5 -lstdc++
GSOURCES=$(wildcard ${DIR_SRC}/*.cu)
ifeq ($(mtally),true)
ifeq ($(dim),1d)
EXENAME=$(DIR_BIN)/box1d_Tmat
else
EXENAME=$(DIR_BIN)/box3d_Tmat
endif
endif

ifeq ($(ctally),true)
ifeq ($(dim),1d)
EXENAME=$(DIR_BIN)/box1d
else
EXENAME=$(DIR_BIN)/box3d
endif
endif

ifeq ($(ftally),true)
ifeq ($(dim),1d)
EXENAME=$(DIR_BIN)/boxf1d
else
EXENAME=$(DIR_BIN)/boxf3d
endif
endif

EXECUTABLE=$(EXENAME)

CSOURCES=$(wildcard ${DIR_SRC}/*.cc)
MSOURCES=$(wildcard ${DIR_SRC}/*.cxx)
COBJECTS=$(patsubst %.cc, ${DIR_OBJ}/%.obj, $(notdir ${CSOURCES}))
MOBJECTS=$(patsubst %.cxx, ${DIR_OBJ}/%.ob, $(notdir ${MSOURCES}))
GOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${GSOURCES}))
LINKJECT=${DIR_OBJ}/dlink.o      
all: $(EXECUTABLE)
$(EXECUTABLE): $(MOBJECTS) $(COBJECTS) $(GOBJECTS) $(LINKJECT)
	$(CC) $(CIDEN) $^ $(LDFLAGS) -o $@
${DIR_OBJ}/%.obj : ${DIR_SRC}/%.cc
	$(CC) $(CIDEN)            $(CMPTYPE) $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC}/%.cxx
	$(NVCC) $(CIDEN)   $(CMPTYPE) $(NCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	$(NVCC) $(CIDEN)   $(CMPTYPE) $(NCFLAGS)  $^ -o $@
$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
	$(NVCC) $(CIDEN) $(LINKLAG) $^ -o $@
clean :  
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.ob  -exec rm -rf {} \;
