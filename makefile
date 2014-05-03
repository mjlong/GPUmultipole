CC=g++
NVCC = nvcc
NCFLAGS=-g -G -dc -arch=sm_21
CCFLAGS=-c -g -I/opt/hdf5/1.8.10-gnu/include
LINKLAG=-arch=sm_21 -dlink
LDFLAGS=-g  -L/opt/hdf5/1.8.10-gnu/lib/ -L/usr/local/cuda/lib64 -lcudart -lhdf5 
CSOURCES=\
CPUComplex.cc\
hdf5IO.cc\
main.cc
GSOURCES=\
CComplex.cu\
Faddeeva.cu\
multipole.cu\
main.cu
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
	$(NVCC) $(NCFLAGS) $^ -o $@
$(LINKJECT) : $(GOBJECTS)
	$(NVCC) $(LINKLAG) $^ -o $@
remove :
	rm -rf *.o  *.obj *~
clean :  
	rm -rf *.o  *.obj *~ $(EXECUTABLE)
