#MITW = 0
#QUICKW_GLOABAL = 11
#QUICKW_TEXTURE = 12
#QUICKW_CONST   = 13
#WHPW = 2
DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
ifeq ($(ver),debug)
NCFLAGS=-g -G -dc -arch=sm_20 -I${DIR_SRC} -I${DIR_SRC}/wfunction #-Xptxas="-v"
CCFLAGS=-c -g -I/home/jlmiao/opt/hdf5/include 
else
NCFLAGS=-dc -arch=sm_20 -I${DIR_SRC} -I${DIR_SRC}/wfunction #-Xptxas="-v"
CCFLAGS=-c -I/home/jlmiao/opt/hdf5/include 
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
  WSOURCES += $(DIR_SRC)/wfunction/Faddeeva.cu 
  WSOURCES += $(DIR_SRC)/wfunction/QuickW.cu
  EXENAME=$(DIR_BIN)/gpumr_quickwc_double
  endif
endif   
#
ifeq ($(FLOAT),1)
  CMPTYPE = -D __CFLOAT
  EXECUTABLE=$(patsubst %_double, %_float, $(EXENAME))
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
	$(CC)             $(CMPTYPE) $(CCFLAGS) $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/%.cu
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
${DIR_OBJ}/%.o : ${DIR_SRC}/wfunction/%.cu
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
#${DIR_OBJ}/$(LINKJECT) : ${DIR_OBJ}/multipole.o
$(LINKJECT) : $(GOBJECTS) $(WOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
	find ${DIR_BIN} -name gpumr_*   -exec rm -rf {} \;
clean :  
	find ${DIR_OBJ} -name *.o   -exec rm -rf {} \;
	find ${DIR_OBJ} -name *.obj -exec rm -rf {} \;
