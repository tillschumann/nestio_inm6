#include "hdf5mpipp.h"
#include "iostream"
#include "mpi.h"
#include "sstream"

void HDF5mpipp::setSize(int index,int i_size)
{
	hsize_t size[1]={i_size};
	H5Dset_extent(dset_ids[index], size);
}

void HDF5mpipp::write(int data[])
{
	hid_t& dset_id = dset_ids[own_id];
	int& n = ns[own_id];

	hsize_t size[1];
    hsize_t offset[1];
    hsize_t dimsext[1] = {buf_size}; 

	/* Extend the dataset */
	n++;
    size[0] = n*buf_size;
	std::cerr << "size=" << size[0] << std::endl;
    //status = H5Dset_extent (dset_id, size);
	std::cerr << "H5Dset_extent complete" << std::endl;

    /* Select a hyperslab in extended portion of dataset  */
    offset[0] = (n-1)*buf_size;
	hid_t filespace = H5Dget_space (dset_id);
    status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL, dimsext, NULL);  

    /* Define memory space */
    hid_t memspace = H5Screate_simple (RANK, dimsext, NULL); 

	/* Define MPI writing property */
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);//H5FD_MPIO_COLLECTIVE);
	
	//std::cerr << "data=" << data[0] << "," << data[1] << std::endl;
	
    /* Write the data to the extended portion of dataset  */
    status = H5Dwrite (dset_id, H5T_NATIVE_INT, memspace, filespace,
                       plist_id, data);
	H5Pclose(plist_id);
	status = H5Sclose (memspace);
    status = H5Sclose (filespace);
}
void HDF5mpipp::newDataSet(std::string datasetname, int size)
{

	std::cout << "Creating new Datasets for " << datasetname << std::endl;

	hsize_t maxdims[1]={H5S_UNLIMITED};
	hsize_t dims[1]={0};
	hsize_t chunk_dims[1]={buf_size};
	/* Create the data space with unlimited dimensions. */
	
    hid_t filespace=H5Screate_simple (RANK, dims, maxdims);
	/* Modify dataset creation properties, i.e. enable chunking  */
	
    hid_t prop=H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk (prop, RANK, chunk_dims);
	
	//hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /* Create a new dataset within the file using chunk 
           creation properties.  */
	hid_t dset_id_other;
	hid_t dset_id;
	
	int clientscount;
	MPI_Comm_size(MPI_COMM_WORLD, &clientscount);
	int len_str = datasetname.size();
	int *sendbuf_len = &len_str;
	int sendcount = 1;
	int *recvbuf_len = new int[clientscount];
	int recvcount=1;
	//
	std::cout << "recvcount=" << recvcount << std::endl;
	MPI_Allgather(sendbuf_len, sendcount, MPI_INT, recvbuf_len, recvcount, MPI_INT, MPI_COMM_WORLD);
	
	//recvcount_len contains length of each node dataset
	
	//count total length
	int total_length=0;
	
	//max_value is used for buffer size
	int max_value=0;
	for (int i=0; i<clientscount; i++) {
		std::cout << "recvbuf_len=" << recvbuf_len[i] << std::endl;
		total_length+=recvbuf_len[i];
		if (max_value<recvbuf_len[i]) max_value = recvbuf_len[i];
	}

	std::cout << "total_length=" << total_length << std::endl;
	std::cout << len_str << std::endl;
	char* sendbuf = new char[max_value];
	strcpy(sendbuf,datasetname.c_str());
	sendcount = max_value;
	char *recvbuf = new char[clientscount*sendcount];
	std::cout << "MPI_Allgather" << std::endl;
	MPI_Allgather(sendbuf, sendcount, MPI_CHAR, recvbuf, max_value, MPI_CHAR, MPI_COMM_WORLD);
	std::cout << "MPI_Allgather done" << std::endl;
	//recvbuf contains names of each node dataset
	std::cout << "clientscount=" << clientscount << std::endl;
	for (int i=0; i<clientscount; i++) {
		char* name = new char[recvbuf_len[i]+1];
		const int start=max_value*i;
		for (int j=0; j<recvbuf_len[i]; j++)
			name[j] = recvbuf[j+start];
		name[recvbuf_len[i]] = NULL;
		dims[0]=size;
		hid_t filespace = H5Screate_simple (RANK, dims, maxdims);

		dset_id=H5Dcreate2 (file, name, H5T_NATIVE_INT, filespace,
                         H5P_DEFAULT, prop, H5P_DEFAULT);
		dset_ids.push_back(dset_id);
		ns.push_back(0);
		status = H5Sclose (filespace);
	}
	H5Pclose(prop);						 
}
HDF5mpipp::HDF5mpipp(std::string filename) : buf_size(2), RANK(1)
{
	hsize_t maxdims[1]={H5S_UNLIMITED};
	hsize_t dims[1]={0};
	hsize_t chunk_dims[1]={buf_size};
	
	hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	//hid_t fapl_id = H5P_DEFAULT;
    /* Create a new file. If file exists its contents will be overwritten. */
    file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	
	MPI_Comm_rank (MPI_COMM_WORLD, &own_id);
	
	H5Pclose(fapl_id);
}

HDF5mpipp::~HDF5mpipp()
{
	/* Close resources */
	for (int i=0; i<dset_ids.size(); i++)
		status = H5Dclose (dset_ids[i]);
    //status = H5Pclose (prop[0]);
    status = H5Fclose (file);
}

void HDF5mpipp::setBufferSize(int s)
{
	buf_size = s;
}

void HDF5mpipp::updateDatasetSizes()
{
	//
}
