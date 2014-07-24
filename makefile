#MITW = 0
#QUICKW_GLOABAL = 11
#QUICKW_TEXTURE = 12
#QUICKW_CONST   = 13
#WHPW = 2
#FOURIEREXPANSION = 3
#QUICKW FOURIER   = 31
DIR_SRC = ./src
DIR_OBJ = ./obj
ifeq ($(compare),1)
DIR_BIN = ./bin/test
endif
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
ifeq ($(ver),debug)
NCFLAGS=-g -G -dc -arch=sm_20 -I${DIR_SRC} -I${DIR_SRC}/wfunction #-Xptxas="-v"
CCFLAGS=-c -g -I/home/jlmiao/opt/hdf5/include 
DIR_BIN = ./bin/debug
else
NCFLAGS=-dc -arch=sm_20 -I${DIR_SRC} -I${DIR_SRC}/wfunction #-Xptxas="-v"
CCFLAGS=-c -I/home/jlmiao/opt/hdf5/include 
DIR_BIN = ./bin/release
endif
LINKLAG=-arch=sm_20 -dlink
LDFLAGS=-L/home/jlmiao/opt/hdf5/lib/ -L/usr/local/cuda-5.5/lib64 -lcudart -lhdf5 
GSOURCES=$(wildcard ${DIR_SRC}/*.cu)
WSOURCES=
# Faddeeva function implementation 
ifeq ($(WFUN),0)
  W_IDEN = -D __MITW
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu
  EXENAME=$(DIR_BIN)/gpumr_mitw_double
else 
  W_IDEN = -D __SAMPLE
  EXENAME=$(DIR_BIN)/gpumr_sample_double
  ifeq ($(WFUN), 3)
  W_IDEN = -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/gpumr_fourierw_double
  endif
  ifeq ($(WFUN), 31)
  W_IDEN = -D __QUICKWF -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/gpumr_quickwf_double
  endif
  ifeq ($(WFUN), 33)
  W_IDEN = -D __INTERPEXP -D __QUICKWF -D __FOURIERW
  WSOURCES += $(DIR_SRC)/wfunction/fourierw.cu
  EXENAME=$(DIR_BIN)/gpumr_quickexpf_double
  endif
  ifeq ($(WFUN),11)
  W_IDEN = -D __QUICKW -D __QUICKWG
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/gpumr_quickwg_double
  endif 
  ifeq ($(WFUN),12)
  W_IDEN = -D __QUICKW -D __QUICKWT
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/gpumr_quickwt_double
  endif
  ifeq ($(WFUN),13)
  W_IDEN = -D __QUICKW -D __QUICKWC
  #FLOAT=1
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/gpumr_quickwc_double
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
CSOURCES=$(wildcard ${DIR_SRC}/*.cc)
COBJECTS=$(patsubst %.cc, ${DIR_OBJ}/%.obj, $(notdir ${CSOURCES}))
GOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o  , $(notdir ${GSOURCES}))
WOBJECTS=$(patsubst %.cu, ${DIR_OBJ}/%.o  , $(notdir ${WSOURCES}))
LINKJECT=${DIR_OBJ}/dlink.o      
all: $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(GOBJECTS) $(WOBJECTS) $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@
${DIR_OBJ}/%.obj : ${DIR_SRC}/%.cc
	@echo $(epoch)
	$(CC)             $(CMPTYPE) $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	@echo $(epoch)
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/wfunction/%.cu
	@echo $(epoch)
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_BIN} -name gpumr_*   -exec rm -rf {} \;
clean :  
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
