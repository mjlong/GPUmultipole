#CPUGPU COPROCESSOR PATTERNS
#ALLCPU = 1
#W()GPU = 2
#XS_GPU = 3
ifeq ($(COP),1)
  COP_IDEN = -D __ALLCPU
endif
ifeq ($(COP),2)
  COP_IDEN = -D __W__GPU
endif
ifeq ($(COP),3)
  COP_IDEN = -D __XS_GPU
endif
#WFUNCTION METHODS
#MITW = 0
#QUICKW_GLOABAL = 11
#QUICKW_TEXTURE = 12
#QUICKW_CONST   = 13
#WHPW = 2
#FOURIEREXPANSION = 3
#QUICKW FOURIER   = 31
#Directories
#TODO: -use-fast-math;
DIR_SRC = ./src
DIR_SRC_PTX = ./src/ptx
DIR_SRC_OPT = ./src/optix
DIR_OBJ = ./obj
DIR_PTX = ./obj/ptx
DIR_HDF5  = /home/jlmiao/opt/hdf5
DIR_CUDA6 = /usr/local/cuda-6.0
DIR_CUDPP = /home/jlmiao/opt/cudpp-2.1
#Include flags
INC_SRC   = -I${DIR_SRC} -I${DIR_SRC}/wfunction
INC_HDF5  = -I${DIR_HDF5}/include
INC_CUDA6 = -I${DIR_CUDA6}/include
INC_CUDPP = -I${DIR_CUDPP}/include
NCINCFLAGS  = $(INC_SRC) $(INC_CUDA6) $(INC_CUDPP) 
CCINCFLAGS  = $(INC_SRC) $(INC_HDF5)
ifeq ($(compare),1)
DIR_BIN = ./bin/test
endif
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
ifeq ($(ver),debug)
NCFLAGS=-g -G -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
NCPLAGS=-ptx -m64 -arch=sm_20 $(NCINCFLAGS)
CCFLAGS=-c -g                 $(CCINCFLAGS) 
DIR_BIN = ./bin/debug
else
NCFLAGS=      -dc -arch=sm_20 $(NCINCFLAGS)  #-Xptxas="-v"
NCPLAGS=-ptx -m64 -arch=sm_20 $(NCINCFLAGS)
CCFLAGS=-c                    $(CCINCFLAGS) 
DIR_BIN = ./bin/release
endif
LINKLAG=   -dlink -arch=sm_20  
LDFLAGS=-L${DIR_HDF5}/lib/ -L${DIR_CUDA6}/lib64 -L${DIR_CUDPP}/lib/ -lcudpp -lcudart -lhdf5 
GSOURCES=$(wildcard ${DIR_SRC}/*.cu)
WSOURCES=
# Faddeeva function implementation 
ifeq ($(WFUN),0)
  W_IDEN = -D __MITW
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu
  EXENAME=$(DIR_BIN)/cgpumr_mitw_double
else 
  W_IDEN = -D __SAMPLE
  EXENAME=$(DIR_BIN)/cgpumr_sample_double
  ifeq ($(WFUN), 3)
  W_IDEN = -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/cgpumr_fourierw_double
  endif
  ifeq ($(WFUN), 31)
  W_IDEN = -D __QUICKWF -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/cgpumr_quickwf_double
  endif
  ifeq ($(WFUN), 33)
  W_IDEN = -D __INTERPEXP -D __QUICKWF -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/cgpumr_quickexpf_double
  endif
  ifeq ($(WFUN),11)
  W_IDEN = -D __QUICKW -D __QUICKWG
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/cgpumr_quickwg_double
  endif 
  ifeq ($(WFUN),12)
  W_IDEN = -D __QUICKW -D __QUICKWT
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/cgpumr_quickwt_double
  endif
  ifeq ($(WFUN),13)
  W_IDEN = -D __QUICKW -D __QUICKWC
  #FLOAT=1
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/cgpumr_quickwc_double
  endif
endif   
#
ifeq ($(compare),1)
  W_IDEN += -D __PROCESS
endif
ifeq ($(compare),2)
  W_IDEN += -D __TRACK
endif
ifeq ($(compare),3)
  W_IDEN += -D __PLOT
endif
#
#
#NOTE: in current philosophy, QuickW tables should construct complex<double> for all cases
#      and fourierw should always give complex<double> as w_val
#TOSUM: never use FLOAT=1
ifeq ($(FLOAT),1)
  CMPTYPE = -D __CFLOAT
  EXECUTABLE=$(patsubst %_double, %_float, $(EXENAME))
  epoch = 'NEVER use FLOAT'
else
  EXECUTABLE=$(EXENAME)
endif
#
ifeq ($(print),track)
RTMETHOD += -D __PRINTTRACK__=$(printid)
endif
CSOURCES=$(wildcard ${DIR_SRC}/*.cc)
MSOURCES=$(wildcard ${DIR_SRC}/*.cxx)
#MSOURCES = main.cxx; 
COBJECTS=$(patsubst %.cc, ${DIR_OBJ}/%.obj, $(notdir ${CSOURCES}))
MOBJECTS=$(patsubst %.cxx, ${DIR_OBJ}/%.ob, $(notdir ${MSOURCES}))
GOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${GSOURCES}))
WOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o, $(notdir ${WSOURCES}))
LINKJECT=${DIR_OBJ}/dlink.o      
all: $(EXECUTABLE)
$(EXECUTABLE): $(MOBJECTS) $(COBJECTS) $(GOBJECTS) $(WOBJECTS) $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@
${DIR_OBJ}/%.obj : ${DIR_SRC}/%.cc
	@echo $(epoch)
	$(CC)             $(COP_IDEN)$(CMPTYPE) $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.ob : ${DIR_SRC}/%.cxx
	@echo $(epoch)
	$(NVCC) $(W_IDEN) $(COP_IDEN) $(RTMETHOD) $(CMPTYPE) $(NCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	@echo $(epoch)
	$(NVCC) $(W_IDEN) $(COP_IDEN)$(RTMETHOD) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/wfunction/%.cu
	@echo $(epoch)
	$(NVCC) $(W_IDEN) $(COP_IDEN)$(CMPTYPE) $(NCFLAGS)  $^ -o $@
$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
clean :  
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.ob  -exec rm -rf {} \;
