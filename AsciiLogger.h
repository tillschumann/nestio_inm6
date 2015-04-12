#include <fstream>
#include "Multimeter.h"
#include "SpikeDetector.h"
#include "abstract_logger.h"

#ifndef ASCIILOGGER_CLASS
#define ASCIILOGGER_CLASS

class AsciiLogger : public ILogger
{
	private:

		std::ofstream* spike_fs_;
		std::ofstream* multi_fs_;

		void print_id_(std::ostream& os, int gid);
		void print_time_(std::ostream& os, const int& timestamp);
		void print_weight_(std::ostream& os, double weight);
		
	public:
		AsciiLogger(std::string sfn, std::string mfn);
		~AsciiLogger();
		//int newDataSet(const std::string, const int, const int);
		//void setSize(int,int);
		//void setBufferSize(int);
		void updateDatasetSizes(const double& t);
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
		void record_spike(int id, int neuron_id, int timestamp);
		void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		void signup_spike(int id, int neuron_id, int expectedSpikeCount);
		void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames);
};

#endif