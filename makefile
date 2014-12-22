#MITW = 0
#QUICKW_GLOABAL = 11
#QUICKW_TEXTURE = 12
#QUICKW_CONST   = 13
#WHPW = 2
#FOURIEREXPANSION = 3
#QUICKW FOURIER   = 31
#Directories
#TODO: -use-fast-math;
DIR_SRC = .
DIR_OBJ = .
DIR_HDF5  = /home/jlmiao/opt/hdf5
DIR_CUDA6 = /usr/local/cuda-6.0
#Include flags
INC_SRC   = -I${DIR_SRC} -I${DIR_SRC}/wfunction
INC_HDF5  = -I${DIR_HDF5}/include
INC_CUDA6 = -I${DIR_CUDA6}/include
NCINCFLAGS  = $(INC_SRC) $(INC_CUDA6) 
CCINCFLAGS  = $(INC_SRC) $(INC_HDF5) 
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
ifeq ($(ver),debug)
NCFLAGS=-g -G -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
CCFLAGS=-c -g                 $(CCINCFLAGS) 
DIR_BIN = .
else
NCFLAGS=      -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
CCFLAGS=-c                    $(CCINCFLAGS) 
DIR_BIN = .
endif
LINKLAG=   -dlink -arch=sm_20  
LDFLAGS=-L${DIR_HDF5}/lib/ -L${DIR_CUDA6}/lib64 -lcudart -lhdf5 
GSOURCES=$(wildcard ${DIR_SRC}/*.cu)
PSOURCES=$(wildcard ${DIR_SRC_PTX}/*.cu)
WSOURCES=
# Faddeeva function implementation 
W_IDEN =
#
#
#
#NOTE: in current philosophy, QuickW tables should construct complex<double> for all cases
#      and fourierw should always give complex<double> as w_val
#TOSUM: never use FLOAT=1
  EXENAME = cpucpu_
  EXECUTABLE=$(EXENAME)
#
CSOURCES=$(wildcard ${DIR_SRC}/*.cc)
MSOURCES=$(wildcard ${DIR_SRC}/*.cpp)
#MSOURCES = main.cxx; 
COBJECTS=$(patsubst %.cc, ${DIR_OBJ}/%.obj, $(notdir ${CSOURCES}))
MOBJECTS=$(patsubst %.cpp, ${DIR_OBJ}/%.ob, $(notdir ${MSOURCES}))
GOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${GSOURCES}))
WOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${WSOURCES}))
#LINKJECT=${DIR_OBJ}/dlink.o      
all:  $(EXECUTABLE)
$(EXECUTABLE): $(MOBJECTS) $(COBJECTS) $(GOBJECTS) $(WOBJECTS) # $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@
${DIR_OBJ}/%.obj : ${DIR_SRC}/%.cc
	$(CC)             $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC}/%.cpp
	$(NVCC) $(W_IDEN) $(NCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	$(NVCC) $(W_IDEN) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/wfunction/%.cu
	$(NVCC) $(W_IDEN) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC_OPT}/%.cpp
	$(NVCC)           $(NCFLAGS) $^ -o $@
#$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
#	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	find ${DIR_OBJ} -name "*.o"   -exec rm -rf {} \;
	find ${DIR_OBJ} -name "*.obj" -exec rm -rf {} \;
	find ${DIR_OBJ} -name "*.ob"  -exec rm -rf {} \;
clean :  
	find ${DIR_OBJ} -name "*.o"   -exec rm -rf {} \;
	find ${DIR_OBJ} -name "*.obj" -exec rm -rf {} \;
	find ${DIR_OBJ} -name "*.ob"  -exec rm -rf {} \;

