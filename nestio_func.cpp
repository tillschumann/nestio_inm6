/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <cmath>
#include <random>
#include "nestio_func.h"

#include <tr1/random>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/*double nestio::rand2(double var, double mean)
{
  std::tr1::uniform_int<int> unif(1, 52);
  return (sqrt(var*rand()/(double)RAND_MAX)*sgn(rand()/(double)RAND_MAX-0.5))+mean;
}

double nestio::rand2(nestio::Distribution &distribution)
{
  return nestio::rand2(distribution.var, distribution.mean);
}*/

int nestio::getThreadHash(int rank, int thread_num) {
    int num_threads = omp_get_max_threads();
    return rank*num_threads+thread_num;
}
int nestio::getThreadHash() {
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    int thread_num = omp_get_thread_num();
    
    return getThreadHash(rank,thread_num);
}

std::ostream& nestio::operator << (std::ostream &o, const nestio::IDistribution &d)
{
  //o << "mean=" << d.getMean() << " var=" << d.getVar();
  d.write_info(o);
  return o;
}

std::ostream& nestio::operator << (std::ostream &o, const nestio::Configuration &c)
{
  o << "logger: " << c.logger << "\n"
    << "bufferSize: " << c.bufferSize << "\n"
    << "numberOfThreads: " << *(c.numberOfThreads) << "\n" 
    << "numberOfProcesses: " << *(c.numberOfProcesses) << "\n"
    << "numberOfSpikeDetectorsPerThread: " << *(c.numberOfSpikeDetectorsPerThread) << "\n"
    << "numberOfMultimetersPerThread: " << *(c.numberOfMultimetersPerThread) << "\n"
    << "samplingIntervalsOfMeter: " << *(c.samplingIntervalsOfMeter) << "\n"
    << "numberOfValuesWrittenByMeter: " << *(c.numberOfValuesWrittenByMeter) << "\n"
    << "deadTimeSpikeDetector: " << *(c.deadTimeSpikeDetector) << "\n"
    << "deadTimeMultimeters: " << *(c.deadTimeMultimeters) << "\n"
    << "deadTimeDeliver: " << *(c.deadTimeDeliver) << "\n";
    
  return o;
}


std::ostream& nestio::operator << (std::ostream &o, const nestio::Loggers &l)
{
  switch (l) {
    case nestio::SIONLIB:
      o << "SIONLIB";
      break;
    case nestio::SIONLIB_BUFFERED:
      o << "SIONLIB_BUFFERED";
      break;
    #ifdef _SIONLIB_COLL
    case nestio::SIONLIB_COLLECTIVE:
      o << "SIONLIB_COLLECTIVE";
      break;
    #endif
    case nestio::HDF5:
      o << "HDF5";
      break;
    case nestio::oHDF5:
      o << "oHDF5";
      break;
    case nestio::oHDF5_BUFFERED:
      o << "oHDF5_BUFFERED";
      break;
    case nestio::oHDF5_COLLECTIVE:
      o << "oHDF5_COLLECTIVE";
      break;
    case nestio::ASCII:
      o << "ASCII";
      break; 
  }
}