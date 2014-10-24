# Edit the following variables as needed
#HDF_INSTALL = /mnt/hdf/packages/hdf5/v189/Linux_2.6/hdf5-1.8.9-linux-static
HDF_INSTALL = /users/schumann/hdf5-1.8.12/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib
#CC          = scorep --user mpic++
CC          = mpic++
#CFLAGS      = -g -O0 -std=c++0x
CFLAGS      = -std=c++0x -DENABLE_TIMING=0 -g -O0 -fopenmp
#LIB         = -lsz -lz -lm
LIB	= -lhdf5 -lrt -lz -lsz

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5_cpp.a

MPISIONLIB =  `~/sionlib/install/sionlib_linux_gomp/bin/sionconfig --libs --mpi -be -64`
MPISIONCFLAGS = `~/sionlib/install/sionlib_linux_gomp/bin/sionconfig --cflags --mpi -be -64` 

SERSIONLIB =  `~/sionlib/install/sionlib_linux_gomp/bin/sionconfig --libs --ser -be -64`
SERSIONCFLAGS = `~/sionlib/install/sionlib_linux_gomp/bin/sionconfig --cflags --ser -be -64` 

all: hdf5mpipp \
 

hdf5mpipp: sionlib_logger.o NESTProxy.o nestio_func.o hdf5mpipp.o AsciiLogger.o SpikeDetector.o Multimeter.o runNESTProxy.o seriestimer.o stopwatch.o scopetimer.o
	$(CC) $(CFLAGS) -o runNESTProxy runNESTProxy.o nestio_func.o hdf5mpipp.o NESTProxy.o sionlib_logger.o AsciiLogger.o SpikeDetector.o Multimeter.o seriestimer.o stopwatch.o scopetimer.o $(INCLUDE) $(LIBSHDF) $(MPISIONLIB) $(LIB)
	
runNESTProxy.o: runNESTProxy.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c runNESTProxy.cpp $(INCLUDE)
	
nestio_func.o: nestio_func.cpp
	$(CC) $(CFLAGS) -c nestio_func.cpp $(INCLUDE)

hdf5mpipp.o: hdf5mpipp.cpp
	$(CC) $(CFLAGS) -c hdf5mpipp.cpp $(INCLUDE)
NESTProxy.o: NESTProxy.cpp NESTProxy.h
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c NESTProxy.cpp $(INCLUDE)
	
SpikeDetector.o: SpikeDetector.cpp SpikeDetector.h
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c SpikeDetector.cpp $(INCLUDE)
Multimeter.o: Multimeter.cpp Multimeter.h
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c Multimeter.cpp $(INCLUDE)
	
sionlib_logger.o: sionlib_logger.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c sionlib_logger.cpp $(INCLUDE)
AsciiLogger.o: AsciiLogger.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c AsciiLogger.cpp $(INCLUDE)
stopwatch.o: timer/stopwatch.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c timer/stopwatch.cpp $(INCLUDE) 
seriestimer.o: timer/seriestimer.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c timer/seriestimer.cpp $(INCLUDE) 
scopetimer.o: timer/scopetimer.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c timer/scopetimer.cpp $(INCLUDE)
	
readSionFile: readSionFile.cpp
	$(CC) $(CFLAGS) $(SERSIONCFLAGS) -o readSionFile readSionFile.cpp $(INCLUDE) $(SERSIONLIB)
	
test_sionlib: test_sionlib.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -o test_sionlib test_sionlib.cpp $(INCLUDE) $(SIONLIB)
	
openmp-hdf5: openmp-hdf5.c
	mpicc -g -O0 -fopenmp -o openmp-hdf5 openmp-hdf5.c $(INCLUDE) $(LIBSHDF) $(LIB)
	
testRand: testRand.o nestio_func.o
	$(CC) $(CFLAGS) -o testRand testRand.o nestio_func.o $(INCLUDE)
	
testRand.o: testRand.cpp
	$(CC) $(CFLAGS) -c testRand.cpp $(INCLUDE)
	
clean: 
	rm -f *.h5 *.o \
 

.SUFFIXES:.o.c
