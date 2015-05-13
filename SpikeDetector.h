#ifndef SPIKEDETECTOR_CLASS
#define SPIKEDETECTOR_CLASS

#include <random>
#include <cmath>
#include <iostream>
#include <vector>
//#include "NESTProxy.h"
#include "nestio_func.h"
#include "abstract_logger.h"

//double rand_value(double var,double mean);

//template < typename L >
class SpikeDetector {
private:
  
  //nestio::Distribution spikes_dist;
  nestio::IDistribution deadTime;
  ILogger* logger;
  bool isSinup;
  
public:
  int spikedetector_id;
  std::vector<int> neuron_ids;
  std::vector<nestio::IDistribution*> spikes_dists;
  SpikeDetector(const int spikedetector_id, ILogger* logger): 
  spikedetector_id(spikedetector_id),
  logger(logger),
  isSinup(false)
  {
    #ifdef _DEBUG_MODE
    std::cout << "Signup SpikeDetector" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tneuron_id=" << spikedetector_id << std::endl;
    #endif
    //std::cout << "\tspikes_mean=" << spikes_dist.mean<< std::endl;
    //logger->signup_spike(neuron_id, 1000, 1);
    //logger->signup_spike(this,1000);
  };
  void connect2Neuron(int id, nestio::IDistribution* dist)
  {
      if (!isSinup) {
	neuron_ids.push_back(id);
	spikes_dists.push_back(dist);
      }
      else {
	std::cerr << "SpikeDetector error: spikedetector has already signed up" << std::endl;
      }
  }
  
  void singup()
  { 
    for (int i=0; i<neuron_ids.size(); i++) {
      //spikes_dists.at(i).init();
      logger->signup_spike(spikedetector_id, neuron_ids.at(i));
    }
    isSinup = true;
  }
  
  void update(double t, int timestamp)
  {
      for (int n=0; n<neuron_ids.size(); n++) {
	//std::cout << "update SpikeDetector " << neuron_ids.at(n) << std::endl;
	int spikes = (int)spikes_dists.at(n)->getValue();
	std::cout << "spikes=" << spikes << std::endl; 
	for (int i=0; i<spikes; i++) {
	  #ifdef _DEBUG_MODE
	    std::cout << "record_spike" << std::endl;
	  #endif
	    logger->record_spike(spikedetector_id, neuron_ids.at(n), timestamp);
	}
      }
  }
};

#endif