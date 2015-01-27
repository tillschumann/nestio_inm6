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
  //virtual void setSize(int,int) {}
  //virtual void setBufferSize(int) {}
  virtual void updateDatasetSizes(const double& t) {}
  virtual void record_multi(Multimeter* multi, int neuron_id, int timestamp, double* v) {}
  virtual void record_spike(SpikeDetector* spike, int neuron_id, int timestamp) {}
  //virtual void signup_spike(int id, int neuron_id, int expectedsize, int buffer_size) {}
  //virtual void signup_multi(int id, int neuron_id, int exactsize, int buffer_size) {}
  
  virtual void signup_spike(SpikeDetector* spike, int neuron_id, int buf) {}
  virtual void signup_multi(Multimeter* multi, int neuron_id, int buf) {}
  
  virtual void createDatasets() {}
};
  
#endif