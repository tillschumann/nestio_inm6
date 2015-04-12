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
  SionBuffer(): buffer(NULL), ptr(0), max_size(0)
  {}
  SionBuffer(int size):buffer(NULL), ptr(0)
  {
    if (size>0) {
      buffer = new char[size];
      max_size=size;
    }
    max_size=0;
  }
  ~SionBuffer()
  {
    if (buffer!=NULL)
      delete[] buffer;
  }
  void extend(int size)
  {
    if (buffer==NULL) {
      buffer = new char[size];
      max_size=size;
    }
    else {
      if (buffer!=NULL)
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
    else {
      std::cout << "SionBuffer: buffer overflow: ptr=" << ptr << " n=" << n << " max_size=" << max_size << std::endl;
      std::cerr << "SionBuffer: buffer overflow: ptr=" << ptr << " n=" << n << " max_size=" << max_size << std::endl;
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
      if (buffer!=NULL) {
	memcpy(new_buffer,buffer,ptr);
	delete[] buffer;
      }
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
		
		struct Parameters_ {
		    nestio::LoggerType loggerType_;
		    double T_;
		    double Tresolution_;
		    
		    bool overwrite_files_;
		    std::string path_;
		    std::string file_extension_;
		    
		    sion_int64 sion_buffer_size_;
		    int logger_buffer_size_;

		    /**
		    * Set default parameter values.
		    * @param Default file name extension, excluding ".".
		    * @param Default value for withtime property
		    * @param Default value for withgid property
		    */
		    Parameters_(const std::string&, const std::string&, nestio::LoggerType, int, sion_int64);

		    //void get(const AsciiLogger2&, DictionaryDatum&) const;  //!< Store current values in dictionary
		    //void set(const AsciiLogger2&, const Buffers_&, const DictionaryDatum&);  //!< Set values from dicitonary
		  };
		  
		  Parameters_ P_;
		
	public:
		Sionlib_logger(): P_(".",".log", nestio::Standard, 100, 100) {};
		Sionlib_logger(const std::string&, const std::string&, int, sion_int64, nestio::LoggerType);
		~Sionlib_logger();
		
		void syncronize(const double& t);
		void initialize(const double T);
		void finalize();
		
		void record_spike(int spikedetector_id, int neuron_id, int timestamp);
		void record_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		
		void srecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void srecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		
		void brecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void brecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		
		void crecord_spike(int spikedetector_id, int neuron_id, int timestamp);
		void crecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		
		void signup_spike(int id, int neuron_id, int expectedSpikeCount);
		void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames);
		
};

extern std::ostream& operator << (std::ostream &o, const Sionlib_logger &l);

#endif
