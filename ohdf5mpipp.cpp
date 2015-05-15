#include "ohdf5mpipp.h"
#include "iostream"
#include "mpi.h"
#include "sstream"
#include <omp.h>


OHDF5mpipp::OHDF5mpipp(std::string filename, int buf_size, nestio::Logger_type logger_type)
: buf_size(buf_size), RANK(2), logger_type(logger_type)
{
	//Init HDF5 file
	MPI_Comm_size(MPI_COMM_WORLD, &clientscount);
	
	hid_t fapl_id;
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /* Create _a new file. If file exists its contents will be overwritten. */
	file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	
	
	MPI_Comm_rank (MPI_COMM_WORLD, &own_id);
	
	int num_threads=omp_get_max_threads();
	if (logger_type == nestio::Buffered || logger_type == nestio::Collective) {
	  
	  buffer_multi = new oHDF5Buffer;
	  buffer_spike = new oHDF5Buffer;
	  //for (int i=0; i<num_threads;i++) {
	  buffer_multi->extend(buf_size*num_threads);
	  buffer_spike->extend(buf_size*num_threads);
	  //}
	}
	
	H5Pclose(fapl_id);
}

OHDF5mpipp::~OHDF5mpipp()
{
  
      std::cout << ".. hdf5 destructor" << std::endl;
      //Update spike datasets
      
      
      if (logger_type == nestio::Buffered || logger_type == nestio::Collective) {
	  
	  delete buffer_multi;
	  delete buffer_spike;
      }
      
	/* Close resources */
	//if (memtype_multi!=NULL)
	//  status = H5Tclose(memtype_multi);
	//if (memtype_spike!=NULL)
	//  status = H5Tclose(memtype_spike);
	
	
	//for (int i=0; i<spike_datasets.size(); i++)
	//	status = H5Dclose (spike_datasets.at(i).dset_id);

	//for (int i=0; i<multi_datasets.size(); i++)
	//	status = H5Dclose (multi_datasets.at(i).dset_id);
}

void OHDF5mpipp::finalize()
{
  H5Pclose(spike_dataset.plist_id);
  H5Pclose(multi_dataset.plist_id);
  
  status = H5Sclose (spike_dataset.filespace);
  //status = H5Sclose (filespace);
  if (spike_datasets.size()>0)
    status = H5Dclose(spike_dataset.dset_id);
  if (multi_datasets.size()>0)
    status = H5Dclose(multi_dataset.dset_id);
  
  
  status = H5Fclose (file);
}

void OHDF5mpipp::setBufferSize(int s)
{
	buf_size = s;
}

int OHDF5mpipp::predictNewSpikeWindowSize(const double& t, HDF5DataSet &dataset) //missing: consider buffer!
{
  if (t>0) {
    int all_new_entries = (int)((double)dataset.entries / (t / (T_-t)));
    int space_needed = all_new_entries-(spike_dataset.window_size-spike_dataset.window_entries)+buf_size/dataset.sizeof_entry;
    if (space_needed>0)
      return space_needed;
    else
      return 0;
  }
  else
    return 200*spike_datasets.size();
}

/**
	extend the dataset sizes
*/
void OHDF5mpipp::syncronize(const double& t)
{
  updateSpikeDataSets(t);
  
  if (logger_type == nestio::Collective) {
    //const int thread_num = omp_get_thread_num();
    #pragma omp single
    {
      buffer_spike->lock();
      if (buffer_spike->getSize()>0) {
	storeContinuousAnalogSignal(spike_dataset, buffer_spike->read(), buffer_spike->getSize()/spike_dataset.sizeof_entry);
	buffer_spike->clear();
      }
      buffer_spike->unlock();
    }
  }
}

void OHDF5mpipp::updateSpikeDataSets(const double& t)
{
#pragma omp single
  {    
    //std::cout << "\tspike_dataset.nodeoffset=" << spike_dataset.nodeoffset << std::endl;
    //std::cout << "\tspike_dataset.window_entries=" << spike_dataset.window_entries << std::endl;
    //std::cout << "\tspike_dataset.window_size=" << spike_dataset.window_size << std::endl;
    //std::cout << "\tspike_dataset.next_window_size=" << spike_dataset.next_window_size << std::endl;
    
    int new_window_size = 0;
    
    int space_needed = predictNewSpikeWindowSize(t, spike_dataset);
    
    if (space_needed>0 && !spike_dataset.extended2next) {
      new_window_size = space_needed;
      //std::cout << "new_window_size" << std::endl;
      spike_dataset.next_window_size = new_window_size;
      //spike_dataset.all_window_size += new_window_size;
    }
    
    int new_window_sizes[clientscount];
    MPI_Allgather(&new_window_size, 1, MPI_INT, new_window_sizes, 1, MPI_INT, MPI_COMM_WORLD);
    
    if (new_window_size>0) {
      spike_dataset.next_nodeoffset=spike_dataset.all_window_size; // old all_window_size
      for (int i=0; i<clientscount; i++) {
	if (i<own_id)
	  spike_dataset.next_nodeoffset += new_window_sizes[i];
      }
    }
    
    bool extendNecessary=false;
    for (int i=0; i<clientscount; i++) {
      extendNecessary = extendNecessary || new_window_sizes[i] > 0;
      spike_dataset.all_window_size += new_window_sizes[i];
    }
    
    //std::cout << "\tspike_dataset.next_nodeoffset=" << spike_dataset.next_nodeoffset << std::endl;
    //std::cout << "\tspike_dataset.all_window_size=" << spike_dataset.all_window_size << std::endl;
    if (extendNecessary) {
      hsize_t size[2]={spike_dataset.all_window_size,1};
      H5Dset_extent(spike_dataset.dset_id, size);
      status = H5Sclose (spike_dataset.filespace);
      spike_dataset.filespace = H5Dget_space (spike_dataset.dset_id);
      spike_dataset.extended2next=true;
      #ifdef _DEBUG_MODE
      std::cout << "H5Dset_extent:" << spike_dataset.all_window_size << std::endl;
      #endif
    }
  }
    
}

void OHDF5mpipp::record_spike(int id, int neuron_id, int timestamp)
{
  const int thread_num = omp_get_thread_num();
    switch (logger_type) {
      case nestio::Standard:
	srecord_spike(id, neuron_id, timestamp);
	break;
      case nestio::Buffered:
	brecord_spike(id, neuron_id, timestamp);
	break;
      case nestio::Collective:
	crecord_spike(id, neuron_id, timestamp);
	break;
    }
}

void OHDF5mpipp::srecord_spike(int id, int neuron_id, int timestamp)
{
  #ifdef _DEBUG_MODE
  std::cout << "recording: \t" << id << "\t" << neuron_id << "\t" << timestamp << std::endl;
  #endif
  char buffer[3*sizeof(int)];
  memcpy(buffer, &id, sizeof(int));
  memcpy(buffer+sizeof(int), &neuron_id, sizeof(int));
  memcpy(buffer+2*sizeof(int), &timestamp, sizeof(int));
  storeContinuousAnalogSignal(spike_dataset,buffer,1);
}	
void OHDF5mpipp::brecord_spike(int id, int neuron_id, int timestamp)
{
  #ifdef _DEBUG_MODE
  std::cout << "brecording: \t" << id << "\t" << neuron_id << "\t" << timestamp << std::endl;
  #endif
  
  //const int thread_num = omp_get_thread_num();
  
  buffer_spike->lock();
  if (!buffer_spike->isEnoughFreeSpace(spike_dataset.sizeof_entry))
  {
    storeContinuousAnalogSignal(spike_dataset, buffer_spike->read(), buffer_spike->getSize()/spike_dataset.sizeof_entry);
    buffer_spike->clear();
  }
  
  *buffer_spike << id << neuron_id << timestamp;
  buffer_spike->unlock();
}

void OHDF5mpipp::crecord_spike(int id, int neuron_id, int timestamp)
{
  #ifdef _DEBUG_MODE
  std::cout << "crecording: \t" << id << "\t" << neuron_id << "\t" << timestamp << std::endl;
  #endif
  
  //const int thread_num = omp_get_thread_num();
 
  buffer_spike->lock();
  buffer_spike->getEnoughFreeSpace(spike_dataset.sizeof_entry);
  
  *buffer_spike << id << neuron_id << timestamp;
  buffer_spike->unlock();
}

void OHDF5mpipp::storeContinuousAnalogSignal(HDF5DataSet &dataset, char* values, int n)
{
    #pragma omp critical
    {
      if (dataset.window_entries<=dataset.window_size) {
	//
      }
      else {
	spike_dataset.extended2next = false;
      	dataset.window_entries=0;
      	dataset.window_size = dataset.next_window_size;
      	if (dataset.next_nodeoffset>0)
      	  dataset.nodeoffset = dataset.next_nodeoffset;
      	else
      	  dataset.nodeoffset+=dataset.all_window_size;
      }
      //std::cout << "n=" << n << std::endl;
      //std::cout << "dataset.nodeoffset=" << dataset.nodeoffset << std::endl;
      //std::cout << "dataset.window_entries=" << dataset.window_entries << std::endl;
      //std::cout << "dataset.window_size=" << dataset.window_size << std::endl;
      //std::cout << "dataset.all_window_size=" << dataset.all_window_size << std::endl;
      
      
      hsize_t offset[2] = {dataset.nodeoffset+dataset.window_entries,0};
      hsize_t dimsext[2] = {n,1};
    
      
      status = H5Sselect_hyperslab (dataset.filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);  

      /* Define MPI writing property */
      //hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    
      /* Define memory space */
      hid_t memspace = H5Screate_simple (RANK, dimsext, NULL);
      
      /* Write the data to the extended portion of dataset  */
      //status = H5Dwrite (dataset.dset_id, dataset.memtype, memspace, dataset.filespace,
      //		      plist_id, values);
      status = H5Dwrite (dataset.dset_id, dataset.memtype, memspace, dataset.filespace,
		      dataset.plist_id, values);
      status = H5Sclose (memspace);
      
      //H5Pclose(plist_id);
  
      dataset.entries+=n;
      dataset.window_entries+=n;
    }
}


void OHDF5mpipp::record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
    char buffer[3*sizeof(int)+multi_dataset.max_numberOfValues*sizeof(double)];
    memcpy(buffer, &id, sizeof(int));
    memcpy(buffer+sizeof(int), &neuron_id, sizeof(int));
    memcpy(buffer+2*sizeof(int), &timestamp, sizeof(int));
    memcpy(buffer+3*sizeof(int), &data[0], data.size()*sizeof(double));
    storeContinuousAnalogSignal(multi_dataset, buffer, 1);  //TODO more values should be possible
}

void OHDF5mpipp::setNodeOffsetAndAllWindowSize(HDF5DataSet &dataset) {
  int window_sizes[clientscount];
  
  MPI_Allgather(&dataset.window_size, 1, MPI_INT, window_sizes, 1, MPI_INT, MPI_COMM_WORLD);
  
  dataset.nodeoffset=0;
  dataset.all_window_size=0;
  for (int i=0; i<clientscount; i++) {
    if (i<own_id)
      dataset.nodeoffset += window_sizes[i];
    dataset.all_window_size += window_sizes[i];
  }
}

void OHDF5mpipp::initialize(const double T)
{
#pragma omp single
  { 
    T_ = T;
    /*
      * Create the compound datatype for memory. 
      */
    
    //spike
    if (spike_datasets.size()>0) {
    spike_dataset.sizeof_entry = 3*sizeof (int);
    spike_dataset.memtype = H5Tcreate (H5T_COMPOUND, spike_dataset.sizeof_entry);
    status = H5Tinsert (spike_dataset.memtype, "id",0, H5T_NATIVE_INT);
    status = H5Tinsert (spike_dataset.memtype, "neuron id", 1*sizeof(int), H5T_NATIVE_INT);
    status = H5Tinsert (spike_dataset.memtype, "timestamp",2*sizeof(int), H5T_NATIVE_INT);
    
    spike_dataset.window_size=predictNewSpikeWindowSize(0, spike_dataset);
    setNodeOffsetAndAllWindowSize(spike_dataset);
    
    registerHDF5DataSet(spike_dataset, "SpikeDetectors"); 
    
    spike_dataset.next_nodeoffset = 0;
    
    spike_dataset.filespace = H5Dget_space (spike_dataset.dset_id);
    }
    
    //multi
    if (multi_datasets.size()>0) {
      //numberOfValues has to be the max numberOfValues of all nodes
      multi_dataset.max_numberOfValues = 0;
      for (int i=0; i<multi_datasets.size(); i++) {
	if (multi_dataset.max_numberOfValues<multi_datasets[i].head.numberOfValues)
	  multi_dataset.max_numberOfValues = multi_datasets[i].head.numberOfValues;
      }
      int send_buf=multi_dataset.max_numberOfValues;
      MPI_Allreduce(&send_buf, &multi_dataset.max_numberOfValues, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      
      multi_dataset.sizeof_entry = 3*sizeof (int)+multi_dataset.max_numberOfValues*sizeof(double);
      multi_dataset.memtype = H5Tcreate (H5T_COMPOUND, multi_dataset.sizeof_entry);
      status = H5Tinsert (multi_dataset.memtype, "id", 0, H5T_NATIVE_INT);
      status = H5Tinsert (multi_dataset.memtype, "neuron id", sizeof(int), H5T_NATIVE_INT);
      status = H5Tinsert (multi_dataset.memtype, "timestamp", 2*sizeof(int), H5T_NATIVE_INT);
      for (int i=0; i<multi_dataset.max_numberOfValues; i++) {
	std::stringstream ss;
	ss << "V" << i;
	status = H5Tinsert (multi_dataset.memtype, ss.str().c_str(),3*sizeof(int)+i*sizeof(double), H5T_NATIVE_DOUBLE);
      }
      
      multi_dataset.window_size=0;
      for (int i=0; i<multi_datasets.size(); i++) {
	multi_dataset.window_size += T_/multi_datasets[i].interval;
      }
      
      //std::cout << "multi window size=" << multi_dataset.window_size  << std::endl;
      setNodeOffsetAndAllWindowSize(multi_dataset);
      
      registerHDF5DataSet(multi_dataset, "Multimeters");
    }
  }
}

void OHDF5mpipp::registerHDF5DataSet(HDF5DataSet& dataset, char* name)
{
  //hsize_t dimsext[2] = {1,1}; 
  //dataset.memspace = H5Screate_simple (RANK, dimsext, NULL);
  
  int chunk_size = buf_size/dataset.sizeof_entry;
  
  std::cout << "chunk_size=" << chunk_size << std::endl;
  std::cout << "dataset.all_window_size=" << dataset.all_window_size << std::endl;
  
  hsize_t maxdims[2]={H5S_UNLIMITED,1};
  hsize_t dims[2]={dataset.all_window_size, 1};
  hsize_t chunk_dims[2]={5*chunk_size,1};			//numberOfValues is to small
  /* Create the data space with unlimited dimensions. */
  
  dataset.plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (logger_type==nestio::Standard || logger_type==nestio::Buffered)
    H5Pset_dxpl_mpio(dataset.plist_id, H5FD_MPIO_INDEPENDENT);
  else
    H5Pset_dxpl_mpio(dataset.plist_id, H5FD_MPIO_COLLECTIVE);
  
  //hid_t filespace=H5Screate_simple (RANK, dims, maxdims);
  dataset.filespace=H5Screate_simple (RANK, dims, maxdims);
  
  /* Modify dataset creation properties, i.e. enable chunking  */
  
  hid_t prop=H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_chunk (prop, RANK, chunk_dims);
  /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */

  hid_t filetype = H5Tcreate (H5T_COMPOUND, 3*8+dataset.max_numberOfValues*8);
  status = H5Tinsert (filetype, "id", 0, H5T_STD_I64BE);
  status = H5Tinsert (filetype, "neuron id", 8, H5T_STD_I64BE);
  status = H5Tinsert (filetype, "timestamp", 16, H5T_STD_I64BE);
  for (int i=0; i<dataset.max_numberOfValues; i++) {
    std::stringstream ss;
    ss << "V" << i;
    status = H5Tinsert (filetype, ss.str().c_str(), 24+i*8, H5T_IEEE_F64BE); //third argument: offset
  }

  /* Create a new dataset within the file using chunk 
      creation properties.  */
  
  std::cout << "H5Dcreate2 name=" << name << " max_numberOfValues=" << dataset.max_numberOfValues << std::endl;

  dataset.dset_id=H5Dcreate2 (file, name, filetype, dataset.filespace,
	    H5P_DEFAULT, prop, H5P_DEFAULT);
  
  status = H5Pclose(prop);
  status = H5Tclose(filetype);
  //status = H5Sclose (filespace);
}

void OHDF5mpipp::signup_spike(int id, int neuron_id)
{
  
  std::stringstream datasetname_ss;
  datasetname_ss << "spike_" << id << "_neuron_" << neuron_id;
  
  oPrivateDataSet ownDataSet;
  ownDataSet.head.id = id;
  ownDataSet.neuron_id = neuron_id;
  //ownDataSet.head.size = -1;
  ownDataSet.head.numberOfValues = 0;
  strcpy(ownDataSet.head.name, datasetname_ss.str().c_str());
  
  #ifdef _DEBUG_MODE
  std::cout << "signup_spike:" << datasetname_ss.str() << std::endl;
  #endif
  
  spike_datasets.push_back(ownDataSet);
}
void OHDF5mpipp::signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames)
{
  std::stringstream datasetname_ss;
  datasetname_ss << "multi_" << id << "_neuron_" << neuron_id;
  
  oPrivateDataSet ownDataSet;
  ownDataSet.head.id = id;
  ownDataSet.neuron_id = neuron_id;
  //ownDataSet.head.size = (int)T_/sampling_interval;      //TODO: simulationTime / sampling_interval
  ownDataSet.head.numberOfValues = valueNames.size();  //TODO
  ownDataSet.interval = sampling_interval;
  strcpy(ownDataSet.head.name, datasetname_ss.str().c_str());
  #ifdef _DEBUG_MODE
  std::cout << "signup_multi:" << datasetname_ss.str() << std::endl;
  #endif
  
  multi_datasets.push_back(ownDataSet);
}

std::ostream& operator << (std::ostream &o, const OHDF5mpipp &l)
{
  o << "OHDF5mpipp";
  return o;
}
