#include "hdf5.h"
#include <cstring>
#include <string>
#include <vector>

#include "Multimeter.h"
#include "SpikeDetector.h"
#include "abstract_logger.h"

#ifndef OHDF5MPIPP_CLASS
#define OHDF5MPIPP_CLASS

#define DataSet_LENGTH 268 //Byte length of struct DataSet

struct HDF5DataSet {
    bool extended2next;
    int window_size;
    int next_nodeoffset;
    int next_window_size;
    int all_window_size;
    int window_entries;
    int max_numberOfValues;
    int sizeof_entry;
    int nodeoffset;
    int entries;
    hid_t memtype;
    hid_t dset_id;
    hid_t filespace;
    hid_t plist_id;
};

struct oDataSet {
    int id;
    int numberOfValues;
    int size;
    char name[256];
};

struct oPrivateDataSet {
    oDataSet head;
    int neuron_id;
    int entries;
    int buffer_size;
    double interval;
    int type;
    hid_t dset_id;
    hid_t dattr_id;
    hid_t memtype;
};

class oHDF5Buffer
{
private:
  omp_lock_t *o_lock;
  int ptr;
  int max_size;
  char* buffer;
public:
  oHDF5Buffer(): buffer(NULL), ptr(0)
  {
    omp_init_lock(o_lock);
  }
  oHDF5Buffer(int size): ptr(0)
  {
    omp_init_lock(o_lock);
    buffer = new char[size];
    max_size=size;
  }
  ~oHDF5Buffer()
  {
    omp_destroy_lock(o_lock);
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
    //std::cout << "n=" << n << " max_size=" << max_size << " ptr=" << ptr << std::endl;
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
  
  void lock() {
    omp_set_lock(o_lock);
  }
  void unlock() {
    omp_unset_lock(o_lock);
  }
  
  
   char* read()
   {
     return buffer;
   }
};

template < typename T >
oHDF5Buffer& operator<<(oHDF5Buffer& buffer, const T v)
{
  //#pragma omp critical
  //{
    buffer.write((const char*)&v, sizeof(T));
  //}
}

class OHDF5mpipp : public ILogger
{
	private:
		HDF5DataSet multi_dataset;
		HDF5DataSet spike_dataset;
		
		hid_t file;
		herr_t       status;
		hsize_t RANK;
		int own_id;
		int buf_size;
		int clientscount;
		std::vector<hid_t> dset_ids;
		std::vector<int> ns;
		std::vector<int> values;
		
		
		nestio::LoggerType logger_type;
		oHDF5Buffer* buffer_multi;
		oHDF5Buffer* buffer_spike;
		
		double T_;
		
		std::vector<int> multi_ns;
		std::vector<int> spike_ns;
		
		std::vector<oPrivateDataSet> multi_datasets;
		std::vector<oPrivateDataSet> spike_datasets;
		
		std::vector<oDataSet> global_multi_datasets;
		std::vector<oDataSet> global_spike_datasets;
		
		std::vector<int> global_shift_multi;
		std::vector<int> global_shift_spike;
		
		std::vector<int> global_number_multi;
		std::vector<int> global_number_spike;
		
		
		int allocated_mem;
		
		int predictNewSpikeWindowSize(const double& t, HDF5DataSet &dataset);
		void updateSpikeDataSets(const double& t);
		void registerHDF5DataSet(HDF5DataSet &dataset, char* name);
		void setNodeOffsetAndAllWindowSize(HDF5DataSet &dataset);
		
	public:
		OHDF5mpipp() {};
		OHDF5mpipp(std::string, int, nestio::LoggerType logger_type);
		~OHDF5mpipp();
		//void write(int* data);
		//void newDataSet(std::string, const int, const int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes(const double& t);
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
		void record_spike(int id, int neuron_id, int timestamp);
		void srecord_spike(int id, int neuron_id, int timestamp);
		void brecord_spike(int id, int neuron_id, int timestamp);
		void crecord_spike(int id, int neuron_id, int timestamp);
		
		void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		void signup_spike(int id, int neuron_id, int expectedSpikeCount);
		void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames, double simulationTime);
		
		//void signup_spike(SpikeDetector* spike, int neuron_id, int buf);
		//void signup_multi(Multimeter* multi, int neuron_id, int buf);
		
		void storeContinuousAnalogSignal(HDF5DataSet &dataset, char* values, int n);
		
		void syncronize(const double& t);
		void initialize(const double T);
		void finalize();
};


extern std::ostream& operator << (std::ostream &o, const OHDF5mpipp &l);

#endif
