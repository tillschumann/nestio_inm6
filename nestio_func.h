/**
/** NESTIO Functions
/**
/**/
#include <math.h>
#include <random>
#include <iostream>
#include <string>
#include <mpi.h>
#include <omp.h>

#include <tr1/random>

#ifndef NESTIO_FUNC
#define NESTIO_FUNC

namespace nestio
{
  class IDistribution
  {
  public:
    virtual ~IDistribution() {}
    virtual void init(double alpha=1.0) {}
    virtual int getIntValue() {}
    virtual double getValue() {}
    virtual void write_info(std::ostream &o) const {}
  };
  class StandardDistribution : public IDistribution
  {
  private:
    std::tr1::ranlux64_base_01 *eng;
    std::tr1::normal_distribution<double> *normal;
  protected:
    double mean;
    double var;
    unsigned seed;
  public:
    StandardDistribution(): mean(0), var(1), eng(NULL)
    {
      normal = new std::tr1::normal_distribution<double>(mean, sqrt(var));
    }
    StandardDistribution(double mean, double var): mean(mean), var(var), eng(NULL)
    {
      normal = new std::tr1::normal_distribution<double>(mean, sqrt(var));
    }
    ~StandardDistribution()
    {
      delete eng;
      delete normal;
    }
    void init(double alpha=1.0)
    {
      int rank;
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      int num_threads = omp_get_max_threads();
      int thread_num = omp_get_thread_num();
      seed = alpha*(rank*num_threads*100 + thread_num);
      eng = new std::tr1::ranlux64_base_01(seed);
    }
    void set(double i_mean, double i_var)
    {
      mean = i_mean;
      var = i_var;
      delete normal;
      normal = new std::tr1::normal_distribution<double>(mean, sqrt(var));
    }
    int getIntValue()
    {
      return (int)getValue();
    }
    double getValue()
    {
      return (*normal)(*eng);
    }
    double getMean() const 
    {
      return mean;
    }
    double getVar() const
    {
      return var;
    }
    void write_info(std::ostream &o) const
    {
      o << "StandardDistribution_" << mean << "_" << var;
    }
    StandardDistribution & operator= (const StandardDistribution & other)
    {
      if (other.eng != NULL)
	eng = new std::tr1::ranlux64_base_01(other.seed);
      else
	eng = NULL;
      seed = other.seed;
      set(other.mean, other.var);
      mean = other.mean;
      var = other.var;
    } 
  };
  
  class PoissonDistribution : public IDistribution
  {
  private:
    std::tr1::ranlux64_base_01 *eng;
    std::tr1::poisson_distribution<int>* dist;
  protected:
    double lambda;
    unsigned seed;
  public:
    PoissonDistribution(): lambda(1), eng(NULL)
    {
      dist = new std::tr1::poisson_distribution<int>(lambda);
    }
    PoissonDistribution(double l): lambda(l), eng(NULL)
    {
      dist = new std::tr1::poisson_distribution<int>(lambda);
    }
    ~PoissonDistribution()
    {
      delete eng;
      delete dist;
    }
    void init(double alpha=1.0)
    {
      int rank;
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      int num_threads = omp_get_max_threads();
      int thread_num = omp_get_thread_num();
      seed = alpha*(rank*num_threads*100 + thread_num);
      eng = new std::tr1::ranlux64_base_01(seed);
    }
    void set(double l)
    {
      lambda = l;
      delete dist;
      dist = new std::tr1::poisson_distribution<int>(lambda);
    }
    int getIntValue()
    {
      return (int)getValue();
    }
    double getValue()
    {
      return (*dist)(*eng);
    }
    void write_info(std::ostream &o) const
    {
      o << "PoissonDistribution_" << lambda;
    }
    PoissonDistribution & operator= (const PoissonDistribution & other)
    {
      if (other.eng != NULL)
	eng = new std::tr1::ranlux64_base_01(other.seed);
      else
	eng = NULL;
      seed = other.seed;
      dist = other.dist;
      set(other.lambda);
      lambda = other.lambda;
    } 
  };
  
  class BinominalDistribution : public IDistribution
  {
  private:
    std::tr1::ranlux64_base_01 *eng;
    std::tr1::binomial_distribution<int>* dist;
  protected:
    double t;
    double p;
    unsigned seed;
  public:
    BinominalDistribution(): t(1), p(0.5), eng(NULL)
    {
      dist = new std::tr1::binomial_distribution<int>(t,p);
    }
    BinominalDistribution(double t, double p): t(t), p(p), eng(NULL)
    {
      dist = new std::tr1::binomial_distribution<int>(t,p);
    }
    ~BinominalDistribution()
    {
      delete eng;
      delete dist;
    }
    void init(double alpha=1.0)
    {
      int rank;
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      int num_threads = omp_get_max_threads();
      int thread_num = omp_get_thread_num();
      seed = alpha*(rank*num_threads*100 + thread_num);
      eng = new std::tr1::ranlux64_base_01(seed);
    }
    void set(double it, double ip)
    {
      t = it;
      p = ip;
      delete dist;
      dist = new std::tr1::binomial_distribution<int>(t,p);
    }
    int getIntValue()
    {
      return (int)getValue();
    }
    double getValue()
    {
      return (*dist)(*eng);
    }
    void write_info(std::ostream &o) const
    {
      o << "BinominalDistribution_" << t << "_" << p;
    }
    BinominalDistribution & operator= (const BinominalDistribution & other)
    {
      if (other.eng != NULL)
	eng = new std::tr1::ranlux64_base_01(other.seed);
      else
	eng = NULL;
      seed = other.seed;
      dist = other.dist;
      set(other.t, other.p);
      t = other.t;
      p = other.p;
    } 
  };
  
  class FixDoubleValue : public IDistribution
  {
  protected:
    double v;
  public:
    FixDoubleValue(): v(0)
    {}
    FixDoubleValue(double value): v(value)
    {}
    void init(double alpha=1.0) {}
    int getIntValue()
    {
      return (int)v;
    }
    double getValue()
    {
      return v;
    }
    void write_info(std::ostream &o) const
    {
      o << "FixDoubleValue_" << v;
    }
    FixDoubleValue & operator= (const FixDoubleValue & other)
    {
      v = other.v;
    } 
  };
  
  class FixIntValue : public IDistribution
  {
  protected:
    int v;
  public:
    FixIntValue(): v(0)
    {}
    FixIntValue(int value): v(value)
    {}
    
    void init(double alpha=1.0) {}
    int getIntValue()
    {
      return v;
    }
    double getValue()
    {
      return (double)v;
    }
    void write_info(std::ostream &o) const
    {
      o << "FixIntValue_" << v;
    }
    FixIntValue & operator= (const FixIntValue & other)
    {
      v = other.v;
    }
  };
  
  
  struct SimSettings {
   double T;
   double Tresolution;
    
  };
  
  enum Loggers {SIONLIB, SIONLIB_BUFFERED, SIONLIB_COLLECTIVE, HDF5, ASCII};
  
  extern std::ostream& operator << (std::ostream &o, const nestio::Loggers &l);
  
  struct Configuration {
    Loggers logger;
    int bufferSize;
    IDistribution *numberOfThreads;
    IDistribution *numberOfProcesses;
    IDistribution *numberOfSpikeDetectorsPerThread;
    IDistribution *numberOfMultimetersPerThread;
    IDistribution *spikesPerDector;
    IDistribution *samplingIntervalsOfMeter;
    IDistribution *numberOfValuesWrittenByMeter;
    IDistribution *deadTimeSpikeDetector;
    IDistribution *deadTimeMultimeters;
    IDistribution *deadTimeDeliver;
  };
  
  extern std::ostream& operator << (std::ostream &o, const nestio::IDistribution &d);
  extern std::ostream& operator << (std::ostream &o, const nestio::Configuration &c);
  
  extern double rand2(double var, double mean);
  extern double rand2(StandardDistribution &distribution);
};

#endif