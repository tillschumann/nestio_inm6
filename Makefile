# Edit the following variables as needed
#HDF_INSTALL = /mnt/hdf/packages/hdf5/v189/Linux_2.6/hdf5-1.8.9-linux-static
HDF_INSTALL = /users/schumann/hdf5-1.8.12/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib
#CC          = scorep --user mpic++
CC          = mpic++
#CFLAGS      = -g -O0 -std=c++0x -D_DEBUG_MODE=1
CFLAGS      = -std=c++0x -DENABLE_TIMING=0 -g -O0 -fopenmp  -DNESTIOPROXY=1 -DHAVE_MPI=1 -DHAVE_SIONLIB=1
#LIB         = -lsz -lz -lm
LIB	= -lhdf5 -lrt -lz -lsz

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5_cpp.a

MPISIONLIB =  `~/usr/sionlib-1.5.5-svn-1690/bin/sionconfig --libs --ompi -be -64`
MPISIONCFLAGS = `~/usr/sionlib-1.5.5-svn-1690/bin/sionconfig --cflags --ompi -be -64` 

SERSIONLIB =  `~/usr/sionlib-1.5.5-svn-1690/bin/sionconfig --libs --ser -be -64`
SERSIONCFLAGS = `~/usr/sionlib-1.5.5-svn-1690/bin/sionconfig --cflags --ser -be -64` 

all: hdf5mpipp \
 

hdf5mpipp: sionlib_logger.o NESTProxy.o nestio_func.o abstract_logger.o hdf5mpipp.o ohdf5mpipp.o AsciiLogger2.o Multimeter.o SpikeDetector.o runNESTProxy.o seriestimer.o stopwatch.o scopetimer.o
	$(CC) $(CFLAGS) -o runNESTProxy runNESTProxy.o abstract_logger.o nestio_func.o hdf5mpipp.o ohdf5mpipp.o NESTProxy.o sionlib_logger.o AsciiLogger2.o Multimeter.o SpikeDetector.o seriestimer.o stopwatch.o scopetimer.o $(INCLUDE) $(LIBSHDF) $(MPISIONLIB) $(LIB)
	
runNESTProxy.o: runNESTProxy.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c runNESTProxy.cpp $(INCLUDE)
	
abstract_logger.o: abstract_logger.cpp abstract_logger.h
	$(CC) $(CFLAGS) -c abstract_logger.cpp $(INCLUDE)
	
nestio_func.o: nestio_func.cpp
	$(CC) $(CFLAGS) -fPIC -c nestio_func.cpp $(INCLUDE)
ohdf5mpipp.o: ohdf5mpipp.cpp
	$(CC) $(CFLAGS) -c ohdf5mpipp.cpp $(INCLUDE)
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
AsciiLogger2.o: AsciiLogger2.cpp
	$(CC) $(CFLAGS) $(MPISIONCFLAGS) -c AsciiLogger2.cpp $(INCLUDE)
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
	
RandomGenerator.so: RandomGenerator.o nestio_func.o
	g++ $(CFLAGS) -shared RandomGenerator.o nestio_func.o -lpython2.6 -o RandomGenerator.so
	
RandomGenerator.o: RandomGenerator.cpp
	g++ $(CFLAGS) -c -fPIC -I/usr/include/python2.6/ RandomGenerator.cpp
clean: 
	rm -f *.h5 *.o \
 

.SUFFIXES:.o.c
