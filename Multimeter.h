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
  double *values;
  ILogger* logger;
  bool isSinup;
  
public:
  
  int multimeter_id;
  double samlpingInterval;
  int numberOfValues;
  std::vector<int> neuron_ids;
  std::vector<std::string> valueNames;
  
  Multimeter(const int multimeter_id, const double interval, nestio::SimSettings &simSettings, const int numberOfValues, ILogger* logger):
  multimeter_id(multimeter_id),
  samlpingInterval(simSettings.Tresolution),
  numberOfValues(numberOfValues),
  logger(logger),
  isSinup(false)
  {
    std::cout << "Signup multi" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tmultimeter_id=" <<multimeter_id << std::endl;
    std::cout << "\tsamlpingInterval=" << samlpingInterval<< std::endl;
    std::cout << "\tnumberOfValues=" << numberOfValues<< std::endl;
    //logger->signup_multi(multimeter_id, (T-dt)/samlpingInterval, 1);
    //logger->signup_multi(this, 1);
    values = new double[numberOfValues];
    for (int i=0; i<numberOfValues; i++) {
      values[i] = multimeter_id+0.1*i;
      std::stringstream ss;
      ss << "V" << i;
      valueNames.push_back(ss.str());
    }
  }
  
  ~Multimeter()
  {
    delete values;
  }
  
  void connect2Neuron(int id)
  {
    std::cout << "connect2Neuron" << std::endl;
    if (!isSinup)
      neuron_ids.push_back(id);
    else
      std::cout << "Multimeter error: multimeter has already signed up" << std::endl;
  }
  
  void singup()
  {
    for (int i=0; i<neuron_ids.size(); i++) {
      logger->signup_multi(this,neuron_ids.at(i),1);
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
      logger->record_multi(multimeter_id, timestamp, values);
     
    }
  }
};

#endif