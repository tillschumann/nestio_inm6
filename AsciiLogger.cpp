#include "AsciiLogger.h"
#include "iostream"
#include "sstream"
#include <fstream> 

void AsciiLogger::print_id_(std::ostream& os, int gid)
{
    os << gid << ' ';
}

void AsciiLogger::print_time_(std::ostream& os, const int& timestampe)
{
  os << timestampe << ' ';
}

void AsciiLogger::print_weight_(std::ostream& os, double weight)
{
    os << weight << ' ';
}



void AsciiLogger::record_spike(SpikeDetector* spike, int neuron_id, int timestamp)
{
  #pragma omp critical
  {
    print_id_(*spike_fs_, spike->spikedetector_id);
    print_id_(*spike_fs_, neuron_id);
    print_time_(*spike_fs_, timestamp);
    *spike_fs_ << '\n';
  }
}

void AsciiLogger::record_multi(Multimeter* multi, int neuron_id, int timestamp, double* v)
{
  #pragma omp critical
  {
    print_id_(*multi_fs_, multi->multimeter_id);
    print_id_(*multi_fs_, neuron_id);
    print_time_(*multi_fs_, timestamp);
    for (int i=0;i<multi->numberOfValues; i++)
      print_weight_(*multi_fs_, v[i]);
    *multi_fs_ << '\n';
  }
}

void AsciiLogger::signup_spike(SpikeDetector* spike, int neuron_id, int buf)
{
  #pragma omp critical
  {
    *spike_fs_ << spike->spikedetector_id << " " << neuron_id << "\n";
  }
}

void AsciiLogger::signup_multi(Multimeter* multi, int neuron_id, int buf)
{
  #pragma omp critical
  {
    *multi_fs_ << multi->multimeter_id << " " << neuron_id << " " << multi->numberOfValues << " ";
    for (int i=0; i<multi->numberOfValues; i++)
      *multi_fs_ << multi->valueNames[i] << " ";
    *multi_fs_ << "\n";
  }
}

AsciiLogger::AsciiLogger(std::string sfn, std::string mfn)
{	
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  
  
  std::stringstream sfn_new;
  sfn_new << sfn << "_" << rank;
  spike_fs_ = new std::ofstream(sfn_new.str(), std::ofstream::out);
  
  std::stringstream mfn_new;
  mfn_new << mfn << "_" << rank;
  multi_fs_ = new std::ofstream(mfn_new.str(), std::ofstream::out);
}



AsciiLogger::~AsciiLogger()
{
    spike_fs_->close();
    multi_fs_->close();
    
    delete spike_fs_;
    delete multi_fs_;
}

/*void AsciiLogger::setBufferSize(int s)
{
    //
}*/

void AsciiLogger::updateDatasetSizes(const double& t)
{
    //
}
