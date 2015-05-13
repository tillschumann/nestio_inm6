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
  virtual ~ILogger() {};
  /*virtual void updateDatasetSizes(const double& t) {}
  virtual void record_multi(Multimeter* multi, int neuron_id, int timestamp, double* v) {}
  virtual void record_spike(SpikeDetector* spike, int neuron_id, int timestamp) {}  
  virtual void signup_spike(SpikeDetector* spike, int neuron_id, int buf) {}
  virtual void signup_multi(Multimeter* multi, int neuron_id, int buf) {}
  virtual void createDatasets() {}*/
  
  
  
  virtual void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data) {}
  virtual void record_spike(int id, int neuron_id, int timestamp) {}
  
  
  
  virtual void signup_spike(int id, int neuron_id) {}
  virtual void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames) {}
  
  virtual void syncronize(const double& t) {}
  virtual void initialize(const double T) {}
  virtual void finalize() {}
};
  
#endif