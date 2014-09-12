#include "sionlib_logger.h"
#include "iostream"
#include "mpi.h"
#include "sstream"
#include <omp.h>

/**
	data is buffered
*/
/*void Sionlib_logger::single_write(double& t, int& v, const int ptr)
{
	/*values.push_back(v);
	if (values.size()>=buf_size) {
		write(&values[0], ptr);
		values.clear();
	}*/
	/*write(t, (char*)&v, sizeof(int), ptr);
}*/

/*void Sionlib_logger::single_write(double& t, double& v, const int ptr)
{
	write(t, (char*)&v, sizeof(double), ptr);
}*/

template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

void Sionlib_logger::record_spike(int neuron_id, double t)
{
    std::cout << "record_spike" << std::endl;
    sion_ensure_free_space(spike_sid, sizeof(int)+sizeof(double));
    fwrite(&neuron_id, sizeof(int), 1, spike_fileptr);
    fwrite(&t, sizeof(double), 1, spike_fileptr);
    std::cout << "record_spike ..end" << std::endl;
    
    //missing: writing weights
}

void Sionlib_logger::record_multi(int neuron_id, int timestamp, double* v)
{
    int multimeter_id=-1;
    int numberOfValues=0;
    for (int i=0; i<header_multi.NodesCount;i++) {
      if (header_multi.nodes[i].owner_id==neuron_id) {
	multimeter_id = header_multi.nodes[i].owner_id;
	numberOfValues = header_multi.nodes[i].numberOfValues;
      }
    }
    
    sion_ensure_free_space(multi_sid, 4*sizeof(int)+(numberOfValues)*sizeof(double));
    fwrite(&multimeter_id, sizeof(int), 1, multi_fileptr);
    fwrite(&neuron_id, sizeof(int), 1, multi_fileptr);
    fwrite(&timestamp, sizeof(int), 1, multi_fileptr);
    fwrite(&numberOfValues, sizeof(int), 1, multi_fileptr);
    fwrite(v, sizeof(double), numberOfValues, multi_fileptr);
}

void Sionlib_logger::signup_spike(int id, int expectedsize, int buf)
{
  //
}
void Sionlib_logger::signup_multi(int id, int size, int buf)
{
  //
}

void Sionlib_logger::signup_spike(SpikeDetector<Sionlib_logger>* spike, int neuron_id, int buf)
{
#pragma omp critical
  {
    header_spike.nodes[header_spike.NodesCount].id=spike->spikedetector_id;
    header_spike.nodes[header_spike.NodesCount].owner_id=neuron_id;
    header_spike.nodes[header_spike.NodesCount].interval=0;
    header_spike.nodes[header_spike.NodesCount].type=2;
    header_spike.NodesCount++;
  }
}

void Sionlib_logger::signup_multi(Multimeter<Sionlib_logger>* multi, int neuron_id, int buf)
{
#pragma omp critical
  {
    header_multi.nodes[header_multi.NodesCount].id=multi->multimeter_id;
    header_multi.nodes[header_multi.NodesCount].owner_id=neuron_id;
    header_multi.nodes[header_multi.NodesCount].interval=multi->samlpingInterval;
    header_multi.nodes[header_multi.NodesCount].numberOfValues=multi->numberOfValues;   //Work around
    for (int i=0; i<multi->numberOfValues; i++) {
      memcpy(header_multi.nodes[header_multi.NodesCount].valueNames[i],multi->valueNames.at(i).c_str(), min(20,(int)multi->valueNames.at(i).size()));
    }
    header_multi.NodesCount++;
  }
}

void Sionlib_logger::createDatasets()
{
#pragma omp single
  {
  std::cout << "createDatasets" << std::endl;
  int numberOfRecords=-1;
  int startOfBody;
  //
  // create Spike createDatasets
  //
  
  fwrite(&header_spike.NodesCount, sizeof(int), 1, spike_fileptr);
  fwrite(&simSettings.Tstart, sizeof(double), 1, spike_fileptr);
  fwrite(&simSettings.T, sizeof(double), 1, spike_fileptr);
  fwrite(&simSettings.Tresolution, sizeof(double), 1, spike_fileptr); // should be Tresolution
  fwrite(&numberOfRecords, sizeof(int), 1, spike_fileptr);
  
  // !!!!! caution
  startOfBody = 3*sizeof(int)+2*sizeof(double)+header_spike.NodesCount*(2*sizeof(int));
  fwrite(&startOfBody, sizeof(int), 1, spike_fileptr);
  
  /*for (int i=0; i<header_spike.NodesCount;i++) {
    fwrite(&header_spike.nodes[i].id, sizeof(int), 1, spike_fileptr);
    fwrite(&header_spike.nodes[i].owner_id, sizeof(int), 1, spike_fileptr);
  }*/
  
  //
  // create Multidatasets
  //
  fwrite(&header_multi.NodesCount, sizeof(int), 1, multi_fileptr);
  fwrite(&simSettings.Tstart, sizeof(double), 1, multi_fileptr);
  fwrite(&simSettings.T, sizeof(double), 1, multi_fileptr);
  fwrite(&simSettings.Tresolution, sizeof(double), 1, multi_fileptr); // should be Tresolution
  fwrite(&numberOfRecords, sizeof(int), 1, multi_fileptr);
  
  // !!!!! caution
  startOfBody = 3*sizeof(int)+2*sizeof(double)+header_multi.NodesCount*(3*sizeof(int)+sizeof(double));
  for (int i=0; i<header_multi.NodesCount;i++)
    startOfBody += header_multi.NodesCount*header_multi.nodes[i].numberOfValues*sizeof(char)*20;
  
  fwrite(&startOfBody, sizeof(int), 1, multi_fileptr);
  for (int i=0; i<header_multi.NodesCount;i++) {
    fwrite(&header_multi.nodes[i].id, sizeof(int), 1, multi_fileptr);
    //fwrite(&header_multi.nodes[i].owner_id, sizeof(int), 1, multi_fileptr);
    fwrite(&header_multi.nodes[i].interval, sizeof(double), 1, multi_fileptr);
    fwrite(&header_multi.nodes[i].numberOfValues, sizeof(int), 1, multi_fileptr);
    for (int j=0; j<header_multi.nodes[i].numberOfValues; j++) {
      //std::cout << "j=" << j << std::endl;
      fwrite(header_multi.nodes[i].valueNames[j], sizeof(char), 20, multi_fileptr);
    }
      
  }
  }
}

Sionlib_logger::Sionlib_logger(std::string filename, int ibuf_size, nestio::SimSettings &simSettings)
:simSettings(simSettings)
{
	//std::cout << "create logger in thread "<< omp_get_thread_num() << std::endl;
        header_multi.NodesCount=0;
	header_spike.NodesCount=0;
  
	/* SION parameters */
	buf_size = ibuf_size;
	int numFiles, own_id, num_procs;
	MPI_Comm lComm;
	//sion_int64 left, bwrote;
	sion_int32 fsblksize;
	char *newfname=NULL;
	/* MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &own_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	/* open parameters */
	numFiles = 1; fsblksize = -1;
	
	std::stringstream spike_ss, multi_ss;
	spike_ss << "spikes_" << filename;
	multi_ss << "multi_" << filename;
	
	char spike_fname[256],multi_fname[256];
	strcpy(spike_fname, spike_ss.str().c_str());
	strcpy(multi_fname, multi_ss.str().c_str());
	
	/* create a new file */
	spike_sid = sion_paropen_mpi(spike_fname, "bw", &numFiles,
	MPI_COMM_WORLD, &lComm,
	&buf_size, &fsblksize,
	&own_id,
	&spike_fileptr, &newfname);
	
	/* create a new file */
	multi_sid = sion_paropen_mpi(multi_fname, "bw", &numFiles,
	MPI_COMM_WORLD, &lComm,
	&buf_size, &fsblksize,
	&own_id,
	&multi_fileptr, &newfname);
}



Sionlib_logger::~Sionlib_logger()
{
	//data in buffer is ignored @TODO

	sion_parclose_mpi(multi_sid);
	sion_parclose_mpi(spike_sid);
}

void Sionlib_logger::setBufferSize(int s)
{
		buf_size = s;
}

//could be used to call sion_ensure_free_space
void Sionlib_logger::updateDatasetSizes()
{
	//
}
