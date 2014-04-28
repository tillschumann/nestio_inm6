#include "sionlib_logger.h"
#include "iostream"
#include "mpi.h"
#include "sstream"

void Sionlib_logger::single_write(int& v)
{
	values.push_back(v);
	if (values.size()>=buf_size) {
		write(&values[0]);
		values.clear();
	}
}

void Sionlib_logger::write(int data[])
{
	std::stringstream ss;
	ss << prefix << "=";
	for (int i=0; i<buf_size;i++)
		ss << " " << data[i];
	ss << " end";
	ss << "\n";
	sion_ensure_free_space(sid, ss.str().size());
	fwrite(ss.str().c_str(), 1, ss.str().size(), fileptr);
}

void Sionlib_logger::newDataSet(std::string datasetname, int size)
{
	prefix=datasetname;
}

Sionlib_logger::Sionlib_logger(std::string filename) : buf_size(2)
{
	/* SION parameters */
	int numFiles, own_id, num_procs;
	MPI_Comm lComm;
	sion_int64 left, bwrote;
	sion_int32 fsblksize;
	char *newfname=NULL;
	/* MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &own_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	/* open parameters */
	numFiles = 1; fsblksize = -1;
	
	char fname[256];
	strcpy(fname, filename.c_str());
	
	/* create a new file */
	sid = sion_paropen_mpi(fname, "bw", &numFiles,
	MPI_COMM_WORLD, &lComm,
	&buf_size, &fsblksize,
	&own_id,
	&fileptr, &newfname);
}

Sionlib_logger::~Sionlib_logger()
{
	sion_parclose_mpi(sid);
}

void Sionlib_logger::setBufferSize(int s)
{
		buf_size = s;
}


void Sionlib_logger::updateDatasetSizes()
{
	//
}