/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <cmath>
#include <random>
#include "nestio_func.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double nestio::rand2(double var, double mean)
{
  return (sqrt(var*rand()/(double)RAND_MAX)*sgn(rand()/(double)RAND_MAX-0.5))+mean;
}

double nestio::rand2(nestio::Distribution &distribution)
{
  return nestio::rand2(distribution.var, distribution.mean);
}

std::ostream& nestio::operator << (std::ostream &o, const nestio::Distribution &d)
{
  o << "mean=" << d.mean << " var=" << d.var;
  return o;
}

std::ostream& nestio::operator << (std::ostream &o, const nestio::Configuration &c)
{
  o << "numberOfThreads: " << c.numberOfThreads << "\n" 
    << "numberOfProcesses: " << c.numberOfProcesses << "\n"
    << "numberOfSpikeDetectorsPerThread: " << c.numberOfSpikeDetectorsPerThread << "\n"
    << "numberOfMultimetersPerThread: " << c.numberOfMultimetersPerThread << "\n"
    << "samlpingIntervalsOfMeter: " << c.samlpingIntervalsOfMeter << "\n"
    << "numberOfValuesWrittenByMeter: " << c.numberOfValuesWrittenByMeter << "\n"
    << "deadTimeSpikeDetector: " << c.deadTimeSpikeDetector << "\n"
    << "deadTimeMultimeters: " << c.deadTimeMultimeters << "\n"
    << "deadTimeDeliver: " << c.deadTimeDeliver << "\n";
    
  return o;
}