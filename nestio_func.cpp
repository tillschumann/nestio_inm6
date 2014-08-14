/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <random>
#include "nestio_func.h"

double nestio::rand2(double var, double mean)
{
  return sqrt(var*(rand()/RAND_MAX-0.5)*2.)+mean;
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
    << "numberOfSpikeDetectorsPerProcess: " << c.numberOfSpikeDetectorsPerProcess << "\n"
    << "numberOfMultimetersPerProcess: " << c.numberOfMultimetersPerProcess << "\n"
    << "samlpingIntervalsOfMeter: " << c.samlpingIntervalsOfMeter << "\n"
    << "numberOfValuesWrittenByMeter: " << c.numberOfValuesWrittenByMeter << "\n"
    << "deadTimeSpikeDetector: " << c.deadTimeSpikeDetector << "\n"
    << "deadTimeMultimeters: " << c.deadTimeMultimeters << "\n";
    
  return o;
}