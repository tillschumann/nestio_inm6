#ifndef ASCIILOGGER_CLASS
#define ASCIILOGGER_CLASS

#include <fstream>
//#include "multimeter.h"
//#include "spike_detector.h"
#include "abstract_logger.h"
#include "recording_device.h"

/**
 * @file AsciiLogger.h
 * Declarations for class AsciiLogger.
 */

namespace nest
{
  class AsciiLogger : public ILogger
  {
	  private:
		  struct sd
		  {
		    int id;
		    int neuron_id;
		    int buf;
		  };
		  struct mm
		  {
		    int id;
		    int neuron_id;
		    int buf;
		  };
		  
		  std::vector<sd> spikedetectors;
		  std::vector<mm> multimeters;

		  std::ofstream spike_fs_;
		  std::ofstream multi_fs_;
		  
		  int n;
		  
		  
	  public:
		  AsciiLogger();
		  ~AsciiLogger();
		  //int newDataSet(const std::string, const int, const int);
		  //void setSize(int,int);
		  //void setBufferSize(int);
		  //void single_write(double& t, int& v, const int ptr);
		  //void single_write(double& t, double& v, const int ptr);
		  
		  void record_spike(int id, int neuron_id, int timestamp);
		  void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		  //void append_value_to_multi_record(int id, int neuron_id, double v, bool endrecord);
		  void signup_spike(int id, int neuron_id, int expectedSpikeCount);
		  void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames);
		  
		  void syncronize(const double t);
		  void initialize(const double T);
		  void finalize();
  };
}

#endif