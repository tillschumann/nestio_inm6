#include <fstream>

#ifndef ASCIILOGGER_CLASS
#define ASCIILOGGER_CLASS

class AsciiLogger
{
	private:

		std::ofstream* spike_fs_;
		std::ofstream* multi_fs_;

		void print_id_(std::ostream& os, int gid);
		void print_time_(std::ostream& os, const double& t);
		void print_weight_(std::ostream& os, double weight);
		
	public:
		AsciiLogger(std::string, int);
		~AsciiLogger();
		//int newDataSet(const std::string, const int, const int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes();
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
		void record_spike(int neuron_id, double t);
		void record_multi(int multimeter_id, double t, double v);
		void signup_spike(int id, int expectedsize);
		void signup_multi(int id, int size);
};

#endif