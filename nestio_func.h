/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <random>
#include <iostream>

#include <tr1/random>

#ifndef NESTIO_FUNC
#define NESTIO_FUNC

namespace nestio
{
  class Distribution
  {
  private:
    std::tr1::ranlux64_base_01 eng;
    std::tr1::normal_distribution<double> normal;
  public:
    Distribution(): mean(0), var(1), normal(0, 1)
    {}
    Distribution(double mean, double var): mean(mean), var(var), normal(mean, sqrt(var))
    {}
    double getValue()
    {
      return normal(eng);
    }
    double mean;
    double var;
    
    Distribution & operator= (const Distribution & other)
    {
      eng = other.eng;
      normal = other.normal;
      mean = other.mean;
      var = other.var;
    } 
  };
  
  
  struct SimSettings {
   double Tstart;
   double T;
   double Tresolution;
    
  };
  
  struct Configuration {
    int numberOfThreads;
    int numberOfProcesses;
    int numberOfSpikeDetectorsPerThread;
    int numberOfMultimetersPerThread;
    Distribution spikesPerDector;
    Distribution samplingIntervalsOfMeter;
    Distribution numberOfValuesWrittenByMeter;
    Distribution deadTimeSpikeDetector;
    Distribution deadTimeMultimeters;
    Distribution deadTimeDeliver;
  };
  
  extern std::ostream& operator << (std::ostream &o, const nestio::Distribution &d);
  extern std::ostream& operator << (std::ostream &o, const nestio::Configuration &c);
  
  extern double rand2(double var, double mean);
  extern double rand2(Distribution &distribution);
};

#endif