#include "hdf5mpipp.h"
#include "iostream"
#include "mpi.h"
#include "sstream"
#include <omp.h>


HDF5mpipp::HDF5mpipp(std::string filename, int buf_size, nestio::SimSettings& simSettings)
: buf_size(buf_size), RANK(2), simSettings(simSettings)
{
	//Init HDF5 file
	MPI_Comm_size(MPI_COMM_WORLD, &clientscount);
	
	hid_t fapl_id;
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /* Create _a new file. If file exists its contents will be overwritten. */
	file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	
	MPI_Comm_rank (MPI_COMM_WORLD, &own_id);
	
	H5Pclose(fapl_id);
}

HDF5mpipp::~HDF5mpipp()
{
  
      std::cout << ".. hdf5 destructor" << std::endl;
      //Update spike datasets
	const int ownNumber = private_ptr_spike_datasets.size();
	const int totalNumber = spike_datasets.size();
	
	if (totalNumber > 0) {
	  int send_newmemsize[ownNumber];
	  for (int i=0;i<ownNumber;i++) {
	    send_newmemsize[i] = spike_datasets[private_ptr_spike_datasets[i]].entries;// * spike_datasets[private_ptr_spike_datasets[i]].buffer_size;
	    std::cout << i << ": send_newmemsize=" << send_newmemsize[i] << std::endl;
	    
	  }
	  int recv_newmemsize[totalNumber];
	  MPI_Allgatherv(send_newmemsize, ownNumber, MPI_INT, recv_newmemsize, &global_number_spike[0], &global_shift_spike[0], MPI_INT, MPI_COMM_WORLD);

	  for (int i=0; i<totalNumber; i++) {
	      std::cout << "new size=" << recv_newmemsize[i] << std::endl;
	      hsize_t size[2]={recv_newmemsize[i],1};
	      H5Dset_extent(spike_datasets.at(i).dset_id, size);
	  }
	}
	/* Close resources */
	for (int i=0; i<private_ptr_multi_datasets.size(); i++) {
	  status = H5Tclose(multi_datasets[private_ptr_multi_datasets[i]].memtype);
	}
	for (int i=0; i<private_ptr_spike_datasets.size(); i++) {
	  status = H5Tclose(spike_datasets[private_ptr_spike_datasets[i]].memtype);
	}
	
	for (int i=0; i<spike_datasets.size(); i++)
		status = H5Dclose (spike_datasets.at(i).dset_id);

	for (int i=0; i<multi_datasets.size(); i++)
		status = H5Dclose (multi_datasets.at(i).dset_id);
    	status = H5Fclose (file);
}

void HDF5mpipp::setBufferSize(int s)
{
	buf_size = s;
}

int HDF5mpipp::predictSpikeMemSpace(const double& t, PrivateDataSet &dataset)
{
    if (dataset.entries+100 > dataset.head.size) {
      //linear interpolation:
      dataset.head.size = dataset.entries+(int)((double)dataset.entries / (t / (simSettings.T-t))) + 100;
      
      //const offset
      //dataset.buffer_size = dataset.entries + 1000;
      return dataset.head.size;
    }
    return -1;
}

/**
	extend the dataset sizes
*/
void HDF5mpipp::updateDatasetSizes(const double& t)
{
	updateSpikeDataSets(t);
}

void HDF5mpipp::updateSpikeDataSets(const double& t)
{
#pragma omp single
  {
   // std::cout << "updateSpikeDataSets" << std::endl;
    const int ownNumber = private_ptr_spike_datasets.size();
    const int totalNumber = spike_datasets.size();
    
    if (totalNumber>0) {
      int send_newmemsize[ownNumber];
      for (int i=0;i<ownNumber;i++)
	send_newmemsize[i] = predictSpikeMemSpace(t,spike_datasets[private_ptr_spike_datasets[i]]);
      
      int recv_newmemsize[totalNumber];
      MPI_Allgatherv(send_newmemsize, ownNumber, MPI_INT, recv_newmemsize, &global_number_spike[0], &global_shift_spike[0], MPI_INT, MPI_COMM_WORLD);

      for (int i=0; i<totalNumber; i++) {
	if (recv_newmemsize[i] > 0) {
	  hsize_t size[2]={recv_newmemsize[i],1};
	  H5Dset_extent(spike_datasets.at(i).dset_id, size);
	}
      }
    }
    //std::cout << "updateSpikeDataSets end" << std::endl;
  }
    
}

void HDF5mpipp::record_spike(SpikeDetector* spike, int neuron_id, int timestamp)
{
   PrivateDataSet* pDataSet;
    //find multidataset
    for (int i=0; i<spike_datasets.size(); i++) {
	if (spike_datasets[i].neuron_id==neuron_id) {
	  pDataSet = &spike_datasets[i];
	}
    }
    //std::cout << neuron_id << ": spike recording: space for " << pDataSet->head.size-pDataSet->entries << " more spikes" << std::endl;
    storeContinuousAnalogSignal(*pDataSet, timestamp, NULL);
    //std::cout << "storeContinuousAnalogSignal end" << std::endl;
}	

void HDF5mpipp::storeContinuousAnalogSignal(PrivateDataSet &pDataSet, int timestamp, double* v)
{
    #pragma omp critical
    {
      hsize_t offset[2] = {pDataSet.entries,0};
      hsize_t dimsext[2] = {1,1}; 
    
      hid_t filespace = H5Dget_space (pDataSet.dset_id);
      status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);  

      /* Define MPI writing property */
      hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      
      char buffer[sizeof(int)+pDataSet.head.numberOfValues*sizeof(double)];
      memcpy(buffer, &timestamp, sizeof(int));
      memcpy(buffer+sizeof(int), v,pDataSet.head.numberOfValues*sizeof(double));
    
      /* Define memory space */
      hid_t memspace = H5Screate_simple (RANK, dimsext, NULL);
      
      /* Write the data to the extended portion of dataset  */
      status = H5Dwrite (pDataSet.dset_id, pDataSet.memtype, memspace, filespace,
		      plist_id, buffer);
      status = H5Sclose (memspace);
      
      H5Pclose(plist_id);
      status = H5Sclose (filespace);
  
      pDataSet.entries++;
    }
}


void HDF5mpipp::record_multi(int neuron_id, int timestamp, double *v)
{
    PrivateDataSet* pDataSet;
    //find multidataset
    for (int i=0; i<multi_datasets.size(); i++) {
	if (multi_datasets[i].neuron_id==neuron_id) {
	  pDataSet = &multi_datasets[i];
	}
    }
    storeContinuousAnalogSignal(*pDataSet, timestamp, v);  //TODO more values should be possible
}

void HDF5mpipp::createDatasets()
{
#pragma omp single
  {
    distributeDatasets(multi_datasets, global_multi_datasets, global_shift_multi, global_number_multi, private_ptr_multi_datasets);
    distributeDatasets(spike_datasets, global_spike_datasets, global_shift_spike, global_number_spike, private_ptr_spike_datasets);

    for (int i=0; i<private_ptr_multi_datasets.size(); i++) {
      
      PrivateDataSet &ownDataSet = multi_datasets[private_ptr_multi_datasets[i]];
      /*
      * Write attributes
      */
      hid_t attr = ownDataSet.dattr_id;
      H5Awrite(attr, H5T_NATIVE_INT, &ownDataSet.neuron_id);
      H5Aclose(attr);    
      
      /*
      * Create the compound datatype for memory. 
      */
      hid_t memtype = H5Tcreate (H5T_COMPOUND, sizeof (int)+ownDataSet.head.numberOfValues*sizeof(double));
      status = H5Tinsert (memtype, "timestamp",0, H5T_NATIVE_INT);
      for (int i=0; i<ownDataSet.head.numberOfValues; i++) {
	std::stringstream ss;
	ss << "V" << i;
	status = H5Tinsert (memtype, ss.str().c_str(),sizeof(int)+i*sizeof(double), H5T_NATIVE_DOUBLE);
      }
      
      ownDataSet.memtype = memtype;
      
    }
    for (int i=0; i<private_ptr_spike_datasets.size(); i++) {
      PrivateDataSet &ownDataSet = spike_datasets[private_ptr_spike_datasets[i]];
      
      hid_t attr = ownDataSet.dattr_id;
      H5Awrite(attr, H5T_NATIVE_INT, &spike_ptrs[i]->neuron_ids.at(0));     //TODO
      H5Aclose(attr);
      
      
      /*
      * Create the compound datatype for memory. 
      */
      hid_t memtype = H5Tcreate (H5T_COMPOUND, sizeof (int)+ownDataSet.head.numberOfValues*sizeof(double));
      status = H5Tinsert (memtype, "timestamp",0, H5T_NATIVE_INT);
      for (int i=0; i<ownDataSet.head.numberOfValues; i++) {
	std::stringstream ss;
	ss << "V" << i;
	status = H5Tinsert (memtype, ss.str().c_str(),sizeof(int)+i*sizeof(double), H5T_NATIVE_DOUBLE);
      }
      
      ownDataSet.memtype = memtype;
    }
  }
}

void HDF5mpipp::distributeDatasets(std::vector<PrivateDataSet> &datasets,std::vector<DataSet> &global_datasets, std::vector<int> &shift, std::vector<int> &global_number, std::vector<int> &private_ptr_datasets)
{
  #ifdef _DEBUG_MODE
  std::cout << "distributeDatasets" << std::endl;
  #endif
  //
  // Gather number of datasets
  //
  int ownNumber = datasets.size();
 
  //std::cout << "ownNumber=" << ownNumber << std::endl;
  int sendcount = 1;
  
  global_number.resize(clientscount);
  int recvcount=1;
  
  MPI_Allgather(&ownNumber, sendcount, MPI_INT, &global_number.at(0), recvcount, MPI_INT, MPI_COMM_WORLD);
  
  //std::cout << "clientscount=" << clientscount << std::endl;
  
  int totalNumber=0;
    for (int i=0;i<clientscount;i++) {
      //std::cout << "global_number[" << i << "]=" << global_number[i] << std::endl;
      totalNumber += global_number[i];
  }
  
  //std::cout << "totalNumber=" << totalNumber << std::endl;
  //
  // Gather all datasetstotalNumberOfMulits
  //
  
  int bytesPerProcess[clientscount];
  int displacement[clientscount];
  displacement[0]=0;
  
  shift.resize(clientscount);
  shift[0]=0;
 
  for (int i=0;i<clientscount;i++) {
    bytesPerProcess[i] = global_number[i]*DataSet_LENGTH;
    //std::cout << "bytesPerProcess[" << i << "]=" << bytesPerProcess[i] << std::endl;
    if (i>0) {
      displacement[i] = displacement[i-1]+bytesPerProcess[i-1];
      shift[i] = displacement[i-1]+global_number[i];
      shift[i] = shift[i-1]+global_number[i];
      //std::cout << "displacement[" << i << "]=" << displacement[i] << std::endl;
    }
  }

  global_datasets.resize(totalNumber);
  
  DataSet send_datasets[ownNumber];
  for (int i=0;i<ownNumber;i++)
    send_datasets[i] = datasets[i].head;
  
  
  //sending structs is not allowed
  // Fix is needed!!!!
  
  MPI_Allgatherv((char*)&send_datasets[0], ownNumber*DataSet_LENGTH, MPI_CHAR, (char*)&global_datasets[0], bytesPerProcess, displacement, MPI_CHAR, MPI_COMM_WORLD);
  
  
  std::vector<PrivateDataSet> datasets_copy = datasets; //workaround
  
  //datasets.resize(totalNumber);
  
  //resize PrivateDataSet list and move old values
  std::vector<PrivateDataSet>::iterator it;
  
  
  PrivateDataSet emptyDataset;
  //std::cout << "own_id=" << own_id << std::endl;
  //std::cout << "clientscount=" << clientscount << std::endl;
  //std::cout << "shift[own_id]=" << shift[own_id] << std::endl;
  //std::cout << "shift[own_id+1]=" << shift[own_id+1] << std::endl;
  //std::cout << "shift[clientscount-1]=" << shift[clientscount-1] << std::endl;
  if (own_id>0) {
    it = datasets.begin();
    datasets.insert(it,shift[own_id],emptyDataset);
  }
  if (own_id+1 < clientscount) {
    it = datasets.end();
    datasets.insert(it,shift[clientscount-1]+global_number[clientscount-1]-shift[own_id+1],emptyDataset);
  }
  
  private_ptr_datasets.resize(ownNumber);
  for (int i=0;i<totalNumber;i++) {
    //create dataset
    bool isPrivateDataset = false;
    for (int j=0;j<ownNumber;j++)
      if (send_datasets[j].id == global_datasets[i].id) {
	isPrivateDataset = true;
	private_ptr_datasets[j] = i;
	//datasets[i] = datasets_copy[j];
      }
    #ifdef _DEBUG_MODE
    std::cout << "id=" <<  global_datasets[i].id << " size=" << global_datasets[i].size << " name=" << global_datasets[i].name << " isPrivateDataset=" << isPrivateDataset<< std::endl;
    #endif
    datasets.at(i).head = global_datasets[i];
    registerHDF5DataSet(datasets.at(i), isPrivateDataset);
  }
}

void HDF5mpipp::registerHDF5DataSet(PrivateDataSet &dataset, bool isPrivateDataset)
{
  
  dataset.entries=0;
  //dataset.buffer_size=buffer_size; //buf_size not used!!!!!!!!!!!!!!!!!!!!
  
  hsize_t maxdims[2]={H5S_UNLIMITED,1};
  hsize_t dims[2]={dataset.head.size, 1};
  hsize_t chunk_dims[2]={buf_size,1};			//numberOfValues is to small
  /* Create the data space with unlimited dimensions. */
  
  hid_t filespace=H5Screate_simple (RANK, dims, maxdims);
  /* Modify dataset creation properties, i.e. enable chunking  */
  
  hid_t prop=H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_chunk (prop, RANK, chunk_dims);  
  /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */

  hid_t filetype = H5Tcreate (H5T_COMPOUND, 8+dataset.head.numberOfValues*8);
  status = H5Tinsert (filetype, "timestamp", 0, H5T_STD_I64BE);
  for (int i=0; i<dataset.head.numberOfValues; i++) {
    std::stringstream ss;
    ss << "V" << i;
    status = H5Tinsert (filetype, ss.str().c_str(), 8+i*8, H5T_IEEE_F64BE); //third argument: offset
  }

  /* Create a new dataset within the file using chunk 
      creation properties.  */

  dataset.dset_id=H5Dcreate2 (file, dataset.head.name, filetype, filespace,
	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // without chunking
  
  
  /* Create attributes */
  #ifdef _DEBUG_MODE
  std::cout << "Create attributes" << std::endl;
  #endif
  
  if (dataset.type == 0) {
    hid_t aid1  = H5Screate(H5S_SCALAR);
    dataset.dattr_id = H5Acreate2(dataset.dset_id, "Neuron ID", H5T_NATIVE_INT, aid1, H5P_DEFAULT, H5P_DEFAULT);
    if (!isPrivateDataset) 
      H5Aclose(dataset.dattr_id);  //only needed if data is written to it later
    H5Sclose(aid1);
  }
  else {
    hid_t aid1  = H5Screate(H5S_SCALAR);
    dataset.dattr_id = H5Acreate2(dataset.dset_id, "Neuron ID", H5T_NATIVE_INT, aid1, H5P_DEFAULT, H5P_DEFAULT);
    if (!isPrivateDataset)
      H5Aclose(dataset.dattr_id); //only needed if data is written to it later
    H5Sclose(aid1);
    
    /*aid1  = H5Screate(H5S_SCALAR);
    dataset.dattr_id = H5Acreate2(dataset.dset_id, "Sampling interval", H5T_NATIVE_INT, aid1, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(aid1);*/
 }
  
  status = H5Tclose(filetype);
  status = H5Sclose (filespace);
}

void HDF5mpipp::signup_spike(SpikeDetector* spike, int neuron_id, int expectedsize)
{  
  #pragma omp critical
  {
    signup_spike(spike->spikedetector_id, neuron_id, expectedsize, 1); //TODO fixed expected size
    
    //not nice!!
    spike_ptrs.push_back(spike);
  }
}

void HDF5mpipp::signup_multi(Multimeter* multi, int neuron_id, int buf)
{
#pragma omp critical
  {
    signup_multi(multi->multimeter_id, neuron_id, (simSettings.T)/multi->samlpingInterval, multi->numberOfValues);
    
    //not nice!!
    multi_ptrs.push_back(multi);
  }
}

void HDF5mpipp::signup_spike(int id, int neuron_id, int expectedsize, int buf)
{
  
  std::stringstream datasetname_ss;
  datasetname_ss << "spike_" << id << "_neuron_" << neuron_id;
  
  PrivateDataSet ownDataSet;
  ownDataSet.head.id = id;
  ownDataSet.neuron_id = neuron_id;
  ownDataSet.head.size = expectedsize;
  ownDataSet.buffer_size = buf;
  ownDataSet.head.numberOfValues = 0;
  ownDataSet.type = 0;
  strcpy(ownDataSet.head.name, datasetname_ss.str().c_str());
  
  #ifdef _DEBUG_MODE
  std::cout << "signup_spike:" << datasetname_ss.str() << std::endl;
  #endif
  
  spike_datasets.push_back(ownDataSet);
}
void HDF5mpipp::signup_multi(int id, int neuron_id, int size, int buf)
{
  std::stringstream datasetname_ss;
  datasetname_ss << "multi_" << id << "_neuron_" << neuron_id;
  
  PrivateDataSet ownDataSet;
  ownDataSet.head.id = id;
  ownDataSet.neuron_id = neuron_id;
  ownDataSet.head.size = size;
  //ownDataSet.buffer_size = buf;
  ownDataSet.head.numberOfValues = buf;
  ownDataSet.type = 1;
  strcpy(ownDataSet.head.name, datasetname_ss.str().c_str());
  #ifdef _DEBUG_MODE
  std::cout << "signup_multi:" << datasetname_ss.str() << std::endl;
  #endif
  
  multi_datasets.push_back(ownDataSet);
}

std::ostream& operator << (std::ostream &o, const HDF5mpipp &l)
{
  o << "HDF5mpipp";
  return o;
}
