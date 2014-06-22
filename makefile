#MITW = 0
#QUCW = 1
#WHPW = 2
CC=h5cc #g++ #h5pcc #g++
NVCC = nvcc
NCFLAGS=-g -G -dc -arch=sm_20
#CCFLAGS=-c -g -I/opt/hdf5/1.8.11-intel/include -I/opt/mpich/3.0.4-intel/include/
CCFLAGS=-c -g -I/home/jlmiao/opt/hdf5/include 
LINKLAG=-arch=sm_20 -dlink
#LDFLAGS=-g -L/opt/mpich/3.0.4-intel/lib/ -L/opt/hdf5/1.8.11-intel/lib/ -L/usr/local/cuda-5.5/lib64 -lcudart -lhdf5 -lmpich
LDFLAGS=-g -L/home/jlmiao/opt/hdf5/lib/ -L/usr/local/cuda-5.5/lib64 -lcudart -lhdf5 
# Faddeeva function implementation 
ifeq ($(WFUN),0)
  W_IDEN = -D __MITW
  GSOURCES=\
  Faddeeva.cu\
  simulation.cu\
  multipole.cu\
  devicemain.cu\
  main.cu
else ifeq ($(WFUN),1)
       W_IDEN = -D __QUICKW
       FLOAT = 1
       GSOURCES=\
       Faddeeva.cu\
       QuickW.cu\
       simulation.cu\
       multipole.cu\
       devicemain.cu\
       main.cu
     else
       W_IDEN = -D __SAMPLE
       GSOURCES=\
       simulation.cu\
       multipole.cu\
       devicemain.cu\
       main.cu
endif   
#
ifeq ($(FLOAT),1)
  CMPTYPE = -D __CFLOAT
endif
#
CSOURCES=\
hdf5IO.cc
COBJECTS=$(CSOURCES:.cc=.obj)
GOBJECTS=$(GSOURCES:.cu=.o)
LINKJECT=dlink.o
EXECUTABLE=tiny
all: $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(GOBJECTS) $(LINKJECT)
	$(CC)  $^ $(LDFLAGS) -o $@

%.obj : %.cc
	$(CC)   $(CCFLAGS) $^ -o $@
%.o : %.cu
	$(NVCC) $(W_IDEN) $(CMPTYPE) $(NCFLAGS)  $^ -o $@
$(LINKJECT) : $(GOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	rm -rf *.o  *.obj *~
clean :  
	rm -rf *.o  *.obj *~ $(EXECUTABLE)
