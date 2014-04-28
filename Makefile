# Edit the following variables as needed
#HDF_INSTALL = /mnt/hdf/packages/hdf5/v189/Linux_2.6/hdf5-1.8.9-linux-static
HDF_INSTALL = /users/schumann/hdf5-1.8.12/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = mpic++
CFLAGS      = -std=c++0x
#LIB         = -lsz -lz -lm
LIB	= -lhdf5 -lrt -lz -lsz

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5_cpp.a

all: hdf5mpipp \
 

hdf5mpipp: NESTProxy.o hdf5mpipp.o runNESTProxy.o
	$(CC) $(CFLAGS) -o runNESTProxy runNESTProxy.o hdf5mpipp.o NESTProxy.o $(INCLUDE) $(LIBSHDF) $(LIB)

runNESTProxy.o: runNESTProxy.cpp
	$(CC) $(CFLAGS) -c runNESTProxy.cpp $(INCLUDE) 

hdf5mpipp.o: hdf5mpipp.cpp
	$(CC) $(CFLAGS) -c hdf5mpipp.cpp $(INCLUDE)
	
NESTProxy.o: NESTProxy.cpp
	$(CC) $(CFLAGS) -c NESTProxy.cpp $(INCLUDE)
clean: 
	rm -f *.h5 *.o \
 

.SUFFIXES:.o.c
