#ifndef MULTIMETER_CLASS
#define MULTIMETER_CLASS

#include <iostream>
#include <vector>
#include <sstream>
#include "nestio_func.h"
#include "abstract_logger.h"

//template < typename L >
class Multimeter {
private:
  double lastRecordT;
  std::vector<double_t> values;
  ILogger* logger;
  bool isSinup;
  
public:
  
  int multimeter_id;
  double samlpingInterval;
  double simulationTime;
  int numberOfValues;
  std::vector<int> neuron_ids;
  std::vector<std::string> valueNames;
  
  Multimeter(const int multimeter_id, const double interval, nestio::SimSettings &simSettings, const int numberOfValues, ILogger* logger):
  multimeter_id(multimeter_id),
  samlpingInterval(simSettings.Tresolution),
  simulationTime(simSettings.T),
  numberOfValues(numberOfValues),
  values(numberOfValues),
  logger(logger),
  lastRecordT(0),
  isSinup(false)
  {
    #ifdef _DEBUG_MODE
    std::cout << "Signup multi" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tmultimeter_id=" <<multimeter_id << std::endl;
    std::cout << "\tsamlpingInterval=" << samlpingInterval<< std::endl;
    std::cout << "\tnumberOfValues=" << numberOfValues<< std::endl;
    #endif
    //logger->signup_multi(multimeter_id, (T-dt)/samlpingInterval, 1);
    //logger->signup_multi(this, 1);
    //values = new double[numberOfValues];
    for (int i=0; i<numberOfValues; i++) {
      values[i] = multimeter_id+0.1*i;
      std::stringstream ss;
      ss << "V" << i;
      valueNames.push_back(ss.str());
    }
  }
  
  ~Multimeter()
  {
  }
  
  void connect2Neuron(int id)
  {
    #ifdef _DEBUG_MODE
    std::cout << "connect2Neuron" << std::endl;
    #endif
    if (!isSinup)
      neuron_ids.push_back(id);
    else
      std::cerr << "Multimeter error: multimeter has already signed up" << std::endl;
  }
  
  void singup()
  {
    for (int i=0; i<neuron_ids.size(); i++) {
      logger->signup_multi(multimeter_id,neuron_ids.at(i),samlpingInterval,valueNames);
    }
    isSinup = true;
  }
  
  void update(double t, int timestamp)
  {
    //std::cout << "update Multimeter " << multimeter_id << std::endl;
    while (t-lastRecordT >= samlpingInterval) {
      lastRecordT+=samlpingInterval;
      //std::cout << "record_multi" << std::endl;
      values[0]+=0.3;
      for (int n=0; n<neuron_ids.size(); n++) {
	values[0]+=0.01;
	logger->record_multi(multimeter_id, neuron_ids.at(n), timestamp, values);
      }
     
    }
  }
};

#endif