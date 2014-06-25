#MITW = 0
#QUICKW_GLOABAL = 11
#QUICKW_TEXTURE = 12
#QUICKW_CONST   = 13
#WHPW = 2
DIR_SRC = .
DIR_OBJ = ./bin
DIR_BIN = ./bin
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
NCFLAGS=-g -G -dc -arch=sm_20 #-Xptxas="-v"
CCFLAGS=-c -g -I/home/jlmiao/opt/hdf5/include 
LINKLAG=-arch=sm_20 -dlink
LDFLAGS=-g -L/home/jlmiao/opt/hdf5/lib/ -L/usr/local/cuda-5.5/lib64 -lcudart -lhdf5 
GSOURCES=\
  main.cu \
  devicemain.cu\
  simulation.cu\
  multipole.cu
# Faddeeva function implementation 
ifeq ($(WFUN),0)
  W_IDEN = -D __MITW
  GSOURCES += Faddeeva.cu
else 
  W_IDEN = -D __SAMPLE
  ifeq ($(WFUN),11)
  W_IDEN = -D __QUICKW -D __QUICKWG
  GSOURCES += Faddeeva.cu QuickW.cu
  endif 
  ifeq ($(WFUN),12)
  W_IDEN = -D __QUICKW -D __QUICKWT
  GSOURCES += Faddeeva.cu QuickW.cu
  endif
  ifeq ($(WFUN),13)
  W_IDEN = -D __QUICKW -D __QUICKWC
  GSOURCES += Faddeeva.cu QuickW.cu
  endif
endif   
#
ifeq ($(FLOAT),1)
  CMPTYPE = -D __CFLOAT
endif
#
CSOURCES= $(wildcard ${DIR_SRC}/*.cc)
COBJECTS=$(CSOURCES:.cc=.obj)
GOBJECTS=$(GSOURCES:.cu=.o)
LINKJECT=dlink.o      
EXECUTABLE=bin/tiny
all: $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(GOBJECTS) $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@

bin/%.obj : %.cc
	$(CC)             $(CMPTYPE) $(CCFLAGS) $< -o $@
bin/%.o : %.cu
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
bin/$(LINKJECT) : $(GOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	rm -rf *.o  *.obj *~
clean :  
	rm -rf *.o  *.obj *~ $(EXECUTABLE)
