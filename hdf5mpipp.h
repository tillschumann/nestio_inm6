#include "hdf5.h"
#include <cstring>
#include <string>
#include <vector>

#include "Multimeter.h"
#include "SpikeDetector.h"
#include "abstract_logger.h"

#ifndef HDF5MPIPP_CLASS
#define HDF5MPIPP_CLASS

#define DataSet_LENGTH 268 //Byte length of struct DataSet
struct DataSet {
    int id;
    int numberOfValues;
    int size;
    char name[256];
};

struct PrivateDataSet {
    DataSet head;
    int neuron_id;
    int entries;
    int buffer_size;
    int type;
    hid_t dset_id;
    hid_t dattr_id;
    hid_t memtype;
};

class HDF5mpipp : public ILogger
{
	private:
		hid_t file;
		herr_t       status;
		hsize_t RANK;
		int own_id;
		int buf_size;
		int clientscount;
		std::vector<hid_t> dset_ids;
		std::vector<int> ns;
		std::vector<int> values;
		
		
		nestio::SimSettings simSettings;
		
		std::vector<int> multi_ns;
		std::vector<int> spike_ns;
		
		std::vector<PrivateDataSet> multi_datasets;
		std::vector<PrivateDataSet> spike_datasets;
		
		std::vector<int> private_ptr_multi_datasets;
		std::vector<int> private_ptr_spike_datasets;
		
		std::vector<DataSet> global_multi_datasets;
		std::vector<DataSet> global_spike_datasets;
		
		std::vector<int> global_shift_multi;
		std::vector<int> global_shift_spike;
		
		std::vector<int> global_number_multi;
		std::vector<int> global_number_spike;
		
		
		std::vector<Multimeter*> multi_ptrs;
		std::vector<SpikeDetector*> spike_ptrs;
		
		
		int allocated_mem;
		
		int predictSpikeMemSpace(const double& t, PrivateDataSet &dataset);
		void updateSpikeDataSets(const double& t);
		void registerHDF5DataSet(PrivateDataSet &dataset, bool isPrivateDataset);
		void distributeDatasets(std::vector<PrivateDataSet> &private_datasets,std::vector<DataSet> &global_datasets, std::vector<int> &shift, std::vector<int> &global_count, std::vector<int> &private_ptr_datasets);
		
	public:
		HDF5mpipp() {};
		HDF5mpipp(std::string, int, nestio::SimSettings&);
		~HDF5mpipp();
		//void write(int* data);
		//void newDataSet(std::string, const int, const int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes(const double& t);
		//void single_write(double& t, int& v, const int ptr);
		//void single_write(double& t, double& v, const int ptr);
		
		void record_spike(int neuron_id, int timestamp);
		void record_multi(int neuron_id, int timestamp, double* v);
		void signup_spike(int id, int neuron_id, int expectedsize, int buffer_size);
		void signup_multi(int id, int neuron_id, int exactsize, int buffer_size);
		
		void signup_spike(SpikeDetector* spike, int neuron_id, int buf);
		void signup_multi(Multimeter* multi, int neuron_id, int buf);
		
		void storeContinuousAnalogSignal(PrivateDataSet &pDataSet, int timestamp, double *v);
		
		void createDatasets();
};


extern std::ostream& operator << (std::ostream &o, const HDF5mpipp &l);

#endif
