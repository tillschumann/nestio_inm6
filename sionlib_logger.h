#include <mpi.h>
#include <omp.h>
#include <sion.h>
#include <vector>
#include "Multimeter.h"
#include "SpikeDetector.h"
#include "abstract_logger.h"

#ifndef SIONLIBLOGGER_CLASS
#define SIONLIBLOGGER_CLASS

struct SionFileHeaderNode
{
  int id;
  int neuron_id;
  int numberOfValues;
  double interval;
  char valueNames[10][20];
};

struct SionFileHeader
{
  int NodesCount;
  std::vector<SionFileHeaderNode> nodes;
  double T;
  double Tresolution;
  int numberOfWrittenData;
};


class SionBuffer
{
private:
  int ptr;
  int max_size;
  char* buffer;
public:
  SionBuffer(): buffer(NULL), ptr(0)
  {}
  SionBuffer(int size): ptr(0)
  {
    buffer = new char[size];
    max_size=size;
  }
  ~SionBuffer()
  {
    delete buffer;
  }
  void extend(int size)
  {
    if (buffer==NULL) {
      buffer = new char[size];
      max_size=size;
    }
    else {
      delete[] buffer;
      buffer = new char[size+max_size];
      max_size+=size;
    }
  }
  void write(const char* v, long unsigned int n)
  {
    if (ptr+n<=max_size) {
      memcpy(buffer+ptr,v,n);
      ptr+=n;
    }
  }
  int getSize()
  {
    return ptr;
  }
  
  void getEnoughFreeSpace(int size)
  {
    if (!isEnoughFreeSpace(size)) {
      int new_max_size = max_size+size*10;
      char* new_buffer = new char[new_max_size];
      max_size=new_max_size;
      memcpy(new_buffer,buffer,ptr);
      delete[] buffer;
      buffer = new_buffer;
    }
  }
  
  bool isEnoughFreeSpace(int size)
  {
    return (ptr+size<max_size);
  }
  
  void clear()
  {
    ptr=0;
  }
  
   char* read()
   {
     return buffer;
   }
};



template < typename T >
SionBuffer& operator<<(SionBuffer& buffer, const T v)
{
  buffer.write((const char*)&v, sizeof(T));
}

class Sionlib_logger : public ILogger
{
	private:
		std::vector<int> spike_sid;
		std::vector<int> multi_sid;
		
		std::vector<SionFileHeader> header_multi;
		std::vector<SionFileHeader> header_spike;
		
		SionBuffer* buffer_multi;
		SionBuffer* buffer_spike;
		
		nestio::SimSettings simSettings;
		
		nestio::LoggerType loggerType;
		
	public:
		Sionlib_logger() {};
		Sionlib_logger(std::string,std::string, sion_int64, nestio::LoggerType, nestio::SimSettings&);
		~Sionlib_logger();
		
		void createDatasets();
		void updateDatasetSizes(const double& t);
		
		void record_spike(SpikeDetector* spike, int neuron_id, int timestamp);
		void record_multi(Multimeter* multi, int neuron_id, int timestamp, double* v);
		
		void srecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void srecord_multi(Multimeter* multi, int neuron_id, int timestamp, double* v);
		
		void brecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void brecord_multi(Multimeter* multi, int neuron_id, int timestamp, double* v);
		
		void crecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void crecord_multi(Multimeter* multi, int neuron_id, int timestamp, double* v);
		
		void signup_spike(SpikeDetector* spike, int neuron_id, int buf);
		void signup_multi(Multimeter* multi, int neuron_id, int buf);
		
};

extern std::ostream& operator << (std::ostream &o, const Sionlib_logger &l);

#endif
