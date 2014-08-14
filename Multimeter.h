#include <iostream>

template < typename L >
class Multimeter {
private:
  int multimeter_id;
  double deadTime_mean;
  double deadTime_var;
  double samlpingInterval;
  double dt;
  double T;
  double lastRecordT;
  L* logger;
  
public:
  
  Multimeter(const int multimeter_id, const double interval, const double dt, const double T, L* logger):
  multimeter_id(multimeter_id),
  samlpingInterval(interval),
  dt(dt),
  T(T),
  logger(logger)
  {
    std::cout << "Signup multi" << std::endl;
    std::cout << "configuration:" << std::endl;
    std::cout << "\tmultimeter_id=" <<multimeter_id << std::endl;
    std::cout << "\tdeadTime_mean=" << deadTime_mean<< std::endl;
    std::cout << "\tdeadTime_var=" << deadTime_var<< std::endl;
    std::cout << "\tsamlpingInterval=" << samlpingInterval<< std::endl;
    std::cout << "\tdt=" <<dt << std::endl;
    logger->signup_multi(multimeter_id, (T-dt)/samlpingInterval, 1);
  }
  
  void update(double t)
  {
    std::cout << "update Multimeter " << multimeter_id << std::endl;
    while (t-lastRecordT >= samlpingInterval) {
      lastRecordT+=samlpingInterval;
      std::cout << "record_multi" << std::endl;
      logger->record_multi(multimeter_id, lastRecordT);
     
    }
  }
};