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
  int type;   //0:int 1:double
  int numberOfValues;
  double interval;
  //std::vector<std::string> valueNames;
  char valueNames[10][20];
};

struct SionFileHeader
{
  int NodesCount;
  SionFileHeaderNode nodes[100];
  double Tstart;
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

namespace nestio
{
  enum SionLoggerType {Standard, Buffered, Collective};
  
};

class Sionlib_logger : public ILogger
{
	private:
		std::vector<int> spike_sid;
		std::vector<int> multi_sid;
		//FILE *spike_fileptr;
		//FILE *multi_fileptr;
		std::vector<int> ids;
		int own_id;
		sion_int64 buf_size;
		std::string prefix;
		std::vector<int> values;
		//std::vecotr<int> datalenghts;
		std::vector<SionFileHeader> header_multi;
		std::vector<SionFileHeader> header_spike;
		
		SionBuffer* buffer_multi;
		SionBuffer* buffer_spike;
		
		nestio::SimSettings simSettings;
		
		nestio::SionLoggerType loggerType;
		
		void write(double& t, char* data, const unsigned int datalength, const int ptr);
		
	public:
		Sionlib_logger() {};
		Sionlib_logger(std::string,std::string, int, nestio::SionLoggerType, nestio::SimSettings&);
		~Sionlib_logger();
		//int newDataSet(const std::string, const int, const int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes(const double& t);
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
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
		
		void createDatasets();
		
};

extern std::ostream& operator << (std::ostream &o, const Sionlib_logger &l);

#endif
