#include "hdf5.h"
#include <cstring>
#include <string>
#include <vector>

#ifndef HDF5MPIPP_CLASS
#define HDF5MPIPP_CLASS

class HDF5mpipp
{
	private:
		hid_t file;
		herr_t       status;
		hsize_t RANK;
		int own_id;
		int buf_size;
		std::vector<hid_t> dset_ids;
		std::vector<int> ns;
		
	public:
		HDF5mpipp(std::string);
		~HDF5mpipp();
		void write(int data[]);
		void newDataSet(std::string, int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes();
};

#endif
