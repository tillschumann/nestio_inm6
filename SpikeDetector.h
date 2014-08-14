#include <random>
#include <cmath>
#include <iostream>
//#include "NESTProxy.h"

//double rand_value(double var,double mean);

template < typename L >
class SpikeDetector {
private:
  int neuron_id;
  double spikes_mean;
  double sqrt_spikes_var;
  double deadTime_mean;
  double deadTime_var;
  L* logger;
  
public:
  
  SpikeDetector(const int neuron_id, const double mean, const double var, L* logger): 
  neuron_id(neuron_id),
  spikes_mean(mean),
  logger(logger)
  {
    std::cout << "Signup SpikeDetector" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tneuron_id=" <<neuron_id << std::endl;
    std::cout << "\tspikes_mean=" << spikes_mean<< std::endl;
    logger->signup_spike(neuron_id, 1000, 1);
    sqrt_spikes_var = sqrt(var);
  };
  
  void update(double t)
  {
      std::cout << "update SpikeDetector " << neuron_id << std::endl;
      int spikes = (int)sqrt(sqrt_spikes_var*rand())+spikes_mean;
      for (int i=0; i<spikes; i++) {
	  std::cout << "record_spike" << std::endl;
	  logger->record_spike(neuron_id, t);
      }
  }
};