CC=g++
H5CC = h5cc
CFLAGS=-c -g -Iinclude 
HFLAGS=-c -g 
LDFLAGS=-lstdc++  
includes = $(wildcard include/*.h)
SOURCES=\
multipole.cpp\
main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
HSOURCES=h5_rdwt.cxx
HOBJECTS=h5_rdwt.o
EXECUTABLE=testmain
all: $(EXECUTABLE)
	@echo "make clean if .h file updated"

$(EXECUTABLE): $(OBJECTS) $(HOBJECTS)
	$(H5CC) $(LDFLAGS) $^ -o $@

%.o : %.cpp ${includes}
	$(CC) $(CFLAGS) $^ -o $@
%.o : %.cxx
	$(CC) $(HFLAGS) $^ -o $@
remove :
	rm -rf *.o  *~
clean :  
	rm -rf *.o  *~ $(EXECUTABLE)
