CC=g++
CFLAGS=-c -g -I/opt/hdf5/1.8.10-gnu/include
LDFLAGS=-L/opt/hdf5/1.8.10-gnu/lib/ -lhdf5 
SOURCES=\
hdf5IO.cc\
main.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=tiny
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC)  $^ $(LDFLAGS) -o $@

%.o : %.cc
	$(CC) $(CFLAGS) $^ -o $@

remove :
	rm -rf *.o  *~
clean :  
	rm -rf *.o  *~ $(EXECUTABLE)
