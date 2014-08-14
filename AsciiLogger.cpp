#include "AsciiLogger.h"
#include "iostream"
#include "sstream"
#include <fstream> 


void AsciiLogger::print_id_(std::ostream& os, int gid)
{
    os << gid << '\t';
}

void AsciiLogger::print_time_(std::ostream& os, const double& t)
{
  os << t << '\t';
}

void AsciiLogger::print_weight_(std::ostream& os, double weight)
{
    os << weight << '\t';
}



void AsciiLogger::record_spike(int neuron_id, double t)
{
    print_id_(*spike_fs_, neuron_id);
    print_time_(*spike_fs_, t);
    *multi_fs_ << 'n';
}

void AsciiLogger::record_multi(int multimeter_id, double t, double v)
{
    print_id_(*multi_fs_, multimeter_id);
    print_time_(*multi_fs_, t);
    print_weight_(*multi_fs_, v);
    *multi_fs_ << 'n';
}

void AsciiLogger::signup_spike(int id, int expectedsize)
{
  //
}

void AsciiLogger::signup_multi(int id, int size)
{
  //
}

AsciiLogger::AsciiLogger(std::string filename, int ibuf_size)
{
	
	std::stringstream spike_ss, multi_ss;
	spike_ss << "spikes_" << filename;
	multi_ss << "multi_" << filename;
	
	spike_fs_ = new std::ofstream(spike_ss.str().c_str(), std::ofstream::out);
	multi_fs_ = new std::ofstream(multi_ss.str(), std::ofstream::out);
}



AsciiLogger::~AsciiLogger()
{
    spike_fs_->close();
    multi_fs_->close();
    
    delete spike_fs_;
    delete multi_fs_;
}

void AsciiLogger::setBufferSize(int s)
{
    //
}

void AsciiLogger::updateDatasetSizes()
{
    //
}
