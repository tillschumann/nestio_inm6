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
#include <vector>

#include <tr1/random>

#ifndef NESTIO_FUNC
#define NESTIO_FUNC


typedef std::string Name;  //TODO: Name class from NEST 
typedef double double_t;   //TODO: double_t class from NEST 
//typedef std::map DictionaryDatum

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
    virtual IDistribution* copy_init() {}
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
      return (int)round(getValue());
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
    IDistribution* copy_init()
    {
	StandardDistribution* copy = new StandardDistribution();
	//std::cout << "this=" << this << std::endl;
	//std::cout << "copy=" << copy << std::endl;
	*copy = *this;
	//std::cout << "this=" << this << std::endl;
	//std::cout << "copy=" << copy << std::endl;
	copy->init();
	return copy;
    }
    
    StandardDistribution & operator= (const StandardDistribution & other)
    {
      //std::cout << "call =" << std::endl;
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
      return (int)round(getValue());
    }
    double getValue()
    {
      return (*dist)(*eng);
    }
    void write_info(std::ostream &o) const
    {
      o << "PoissonDistribution_" << lambda;
    }
    IDistribution* copy_init()
    {
	PoissonDistribution* copy = new PoissonDistribution();
	*copy = *this;
	copy->init();
	return copy;
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
      return (int)round(getValue());
    }
    double getValue()
    {
      return (*dist)(*eng);
    }
    void write_info(std::ostream &o) const
    {
      o << "BinominalDistribution_" << t << "_" << p;
    }
    IDistribution* copy_init()
    {
	BinominalDistribution* copy = new BinominalDistribution();
	*copy = *this;
	copy->init();
	return copy;
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
      return (int)round(v);
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
    IDistribution* copy_init()
    {
	FixDoubleValue* copy = new FixDoubleValue();
	*copy = *this;
	copy->init();
	return copy;
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
    IDistribution* copy_init()
    {
	FixIntValue* copy = new FixIntValue();
	*copy = *this;
	copy->init();
	return copy;
    }
  };
  
  
  struct SimSettings {
   double T;
   double Tresolution;
    
  };
  
  enum Loggers{
    SIONLIB,
    SIONLIB_BUFFERED,
    #ifdef _SIONLIB_COLL
    SIONLIB_COLLECTIVE,
    #endif
    HDF5,
    oHDF5,
    oHDF5_BUFFERED,
    oHDF5_COLLECTIVE,
    ASCII
  };
  enum Logger_type {Standard, Buffered, Collective};
  
  extern std::ostream& operator << (std::ostream &o, const nestio::Loggers &l);
  
  struct Multimeter_Config {
    IDistribution *numberOfValuesWritten;
    IDistribution *samplingIntervals;
    //IDistribution *deadTime;
  };
  
  struct Spikedetector_Config {
    IDistribution *spikesPerDector;
  };
  
  struct Configuration {
    Loggers logger;
    int bufferSize;
    std::map<int,std::vector<Multimeter_Config>> multimeter_configs;
    std::map<int,std::vector<Spikedetector_Config>> spikedetector_configs;
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
  
  extern int getThreadHash();
  extern int getThreadHash(int rank, int thread_num);
  
  /*template <typename T>
  void updateValue(const DictionaryDatum& d, const Names& n, T& v) 
  {
    v = d[n];
  }*/
};

#endif