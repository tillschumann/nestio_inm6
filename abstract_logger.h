#include "Multimeter.h"
#include "SpikeDetector.h"

#ifndef ABSTRACT_LOGGER_CLASS
#define ABSTRACT_LOGGER_CLASS

class SpikeDetector;
class Multimeter;

class ILogger
{
public:
  ILogger() {};
  virtual ~ILogger() {
     #ifdef _DEBUG_MODE
    //std::cout << "Warning: ILogger virtual destructor called" << std::endl;
    #endif
  };
  
  
  virtual void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual record_multi called" << std::endl;
    #endif
  }
  virtual void record_spike(int id, int neuron_id, int timestamp) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual record_spike called" << std::endl;
    #endif
  }
  virtual void signup_spike(int id, int neuron_id) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual signup_spike called" << std::endl;
    #endif
  }
  virtual void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual signup_multi called" << std::endl;
    #endif
  }
  virtual void synchronize(const double& t) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual synchronize called" << std::endl;
    #endif
  }
  virtual void initialize(const double T) {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual initialize called" << std::endl;
    #endif
  }
  virtual void finalize() {
    #ifdef _DEBUG_MODE
    std::cout << "Warning: ILogger virtual finalize called" << std::endl;
    #endif
  }
};
  
#endif