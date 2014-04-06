CC=g++
H5CC = h5cc
CFLAGS=-c -g -Iinclude 
HFLAGS=-c -lstdc++ 
LDFLAGS=  
includes = $(wildcard include/*.h)
SOURCES=\
#CPUComplex.cpp\
#isotope.cpp\
#multipole.cpp\
testmain.cpp
OBJECTS=$(SOURCES:.cpp=.o)
HSOURCES=h5_rdwt.cpp
HOBJECTS=h5_rdwt.o
EXECUTABLE=testmain
all: $(EXECUTABLE)
	@echo "make clean if .h file updated"

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

%.o : %.cpp ${includes}
	$(CC) $(CFLAGS) $^ -o $@
%.o : %.cxx
	$(H5CC) $(HFLAGS) $^ -o $@
remove :
	rm -rf *.o  *~
clean :  
	rm -rf *.o  *~ $(EXECUTABLE)
