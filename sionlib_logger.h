#include <mpi.h>
#include <omp.h>
#include <sion.h>
#include <vector>
#include "Multimeter.h"
#include "SpikeDetector.h"

#ifndef SIONLIBLOGGER_CLASS
#define SIONLIBLOGGER_CLASS

struct SionFileHeaderNode
{
  int id;
  int owner_id;
  int type;   //0:int 1:double
  int numberOfValues;
  double interval;
  //std::vector<std::string> valueNames;
  char valueNames[10][20];
};

struct SionFileHeader
{
  int NodesCount;
  SionFileHeaderNode nodes[10];
  double Tstart;
  double T;
  double Tresolution;
};


class SionBuffer
{
private:
  int dim0;
  int dim1;
  int ptr;
  int max_size;
  char* buffer;
public:
  SionBuffer(int dim0, int dim1): buffer(NULL), ptr(0)
  {}
  ~SionBuffer()
  {
    delete buffer;
  }
  void extend(int n)
  {
    if (buffer==NULL) {
      buffer = new char[n];
      max_size=n;
    }
    else {
      delete buffer;
      buffer = new char[n+max_size];
      max_size+=n;
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
  
  bool isFull()
  {
    return (ptr==max_size);
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

class Sionlib_logger
{
	private:
		int spike_sid;
		int multi_sid;
		FILE *spike_fileptr;
		FILE *multi_fileptr;
		std::vector<int> ids;
		int own_id;
		sion_int64 buf_size;
		std::string prefix;
		std::vector<int> values;
		//std::vecotr<int> datalenghts;
		SionFileHeader header_multi;
		SionFileHeader header_spike;
		
		SionBuffer buffer_multi;
		SionBuffer buffer_spike;
		
		nestio::SimSettings simSettings;
		
		void write(double& t, char* data, const unsigned int datalength, const int ptr);
		
	public:
		Sionlib_logger(std::string, int, nestio::SimSettings&);
		~Sionlib_logger();
		//int newDataSet(const std::string, const int, const int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes();
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
		void record_spike(int neuron_id, double t);
		void record_multi(int multimeter_id, int timestamp, double* v);
		
		void brecord_spike(int neuron_id, double t);
		void brecord_multi(int neuron_id, int timestamp, double* v);
		
		void signup_spike(int id, int expectedsize, int buf);
		void signup_spike(SpikeDetector<Sionlib_logger>* spike, int neuron_id, int buf);
		void signup_multi(int id, int size, int buf);
		void signup_multi(Multimeter<Sionlib_logger>* multi, int neuron_id, int buf);
		
		void createDatasets();
		
};

#endif
