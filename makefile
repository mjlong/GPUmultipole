CC=g++
H5CC = h5cc
CFLAGS=-c -g 
CFLAGS +=-I/opt/hdf5/1.8.10-gnu/include/
CFLAGS +=-I/opt/hdf5/hdf5-1.8.10/src/
HFLAGS=-c -g -lstdc++ 
LDFLAGS=-lstdc++  
includes = 
SOURCES=\
CPUComplex.cpp\
Faddeeva.cpp\
multipole.cpp\
main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
HSOURCES=h5_rdwt.cxx
HOBJECTS=$(HSOURCES:.cxx=.o)  
EXECUTABLE=testmain
all: $(EXECUTABLE)
	@echo "make clean if .h file updated"

$(EXECUTABLE): $(OBJECTS) $(HOBJECTS)
	$(H5CC) $(LDFLAGS) $^ -o $@

%.o : %.cpp 
	$(CC) $(CFLAGS) $^ -o $@
%.o : %.cxx
	$(CC) $(CFLAGS) $^ -o $@
remove :
	rm -rf *.o  *~
clean :  
	rm -rf *.o  *~ $(EXECUTABLE)
