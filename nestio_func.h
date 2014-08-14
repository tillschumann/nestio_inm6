/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <random>
#include <iostream>

#ifndef NESTIO_FUNC
#define NESTIO_FUNC

namespace nestio
{
  struct Distribution {
    double mean;
    double var;
  };
  
  struct Configuration {
    int numberOfThreads;
    int numberOfProcesses;
    int numberOfSpikeDetectorsPerProcess;
    int numberOfMultimetersPerProcess;
    Distribution spikesPerDector;
    Distribution samlpingIntervalsOfMeter;
    Distribution numberOfValuesWrittenByMeter;
    Distribution deadTimeSpikeDetector;
    Distribution deadTimeMultimeters;
  };
  
  extern std::ostream& operator << (std::ostream &o, const nestio::Distribution &d);
  extern std::ostream& operator << (std::ostream &o, const nestio::Configuration &c);
  
  extern double rand2(double var, double mean);
  extern double rand2(Distribution &distribution);
};

#endif