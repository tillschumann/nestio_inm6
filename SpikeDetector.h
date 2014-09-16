#include <random>
#include <cmath>
#include <iostream>
#include <vector>
//#include "NESTProxy.h"
#include "nestio_func.h"

#ifndef SPIKEDETECTOR_CLASS
#define SPIKEDETECTOR_CLASS

//double rand_value(double var,double mean);

template < typename L >
class SpikeDetector {
private:
  
  //nestio::Distribution spikes_dist;
  nestio::Distribution deadTime;
  L* logger;
  bool isSinup;
  
public:
  int spikedetector_id;
  std::vector<int> neuron_ids;
  std::vector<nestio::Distribution> spikes_dists;
  SpikeDetector(const int spikedetector_id, L* logger): 
  spikedetector_id(spikedetector_id),
  logger(logger),
  isSinup(false)
  {
    std::cout << "Signup SpikeDetector" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tneuron_id=" << spikedetector_id << std::endl;
    //std::cout << "\tspikes_mean=" << spikes_dist.mean<< std::endl;
    //logger->signup_spike(neuron_id, 1000, 1);
    //logger->signup_spike(this,1000);
  };
  void connect2Neuron(int id, nestio::Distribution dist)
  {
      if (!isSinup) {
	neuron_ids.push_back(id);
	spikes_dists.push_back(dist);
      }
      else {
	std::cout << "SpikeDetector error: spikedetector has already signed nup" << std::endl;
      }
  }
  
  void singup()
  { 
    for (int i=0; i<neuron_ids.size(); i++)
      logger->signup_spike(this, neuron_ids.at(i),1000);
    isSinup = true;
  }
  
  void update(double t, int timestamp)
  {
      for (int n=0; n<neuron_ids.size(); n++) {
	//std::cout << "update SpikeDetector " << neuron_ids.at(n) << std::endl;
	int spikes = (int)nestio::rand2(spikes_dists.at(n));
	for (int i=0; i<spikes; i++) {
	    //std::cout << "record_spike" << std::endl;
	    logger->record_spike(neuron_ids.at(n), timestamp);
	}
      }
  }
};

#endif