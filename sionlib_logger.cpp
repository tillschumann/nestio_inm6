#include "sionlib_logger.h"
#include "iostream"
#include "mpi.h"
#include "sstream"
#include <omp.h>


template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

void Sionlib_logger::crecord_spike(int neuron_id, double t)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif
 
  buffer_spike[thread_num].getEnoughFreeSpace(sizeof(int)+sizeof(double));
  buffer_spike[thread_num] << neuron_id << t;
}

void Sionlib_logger::brecord_spike(int neuron_id, double t)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif
  
  if (!buffer_spike[thread_num].isEnoughFreeSpace(sizeof(int)+sizeof(double)))
  {
    sion_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
    buffer_spike[thread_num].clear();
  }
  buffer_spike[thread_num] << neuron_id << t;
}

void Sionlib_logger::crecord_multi(int neuron_id, int timestamp, double* v)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif
  
  int multimeter_id=-1;
    int numberOfValues=0;
    for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
      if (header_multi[thread_num].nodes[i].owner_id==neuron_id) {
	multimeter_id = header_multi[thread_num].nodes[i].owner_id;
	numberOfValues = header_multi[thread_num].nodes[i].numberOfValues;
      }
    }
  buffer_multi[thread_num].getEnoughFreeSpace(4*sizeof(int)+numberOfValues*sizeof(double));
  buffer_multi[thread_num] << multimeter_id << neuron_id << timestamp << numberOfValues;
  for (int i=0; i<numberOfValues; i++)
  {
    buffer_multi[thread_num] << v[i];
  }
}

void Sionlib_logger::brecord_multi(int neuron_id, int timestamp, double* v)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif
  
  int multimeter_id=-1;
    int numberOfValues=0;
    for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
      if (header_multi[thread_num].nodes[i].owner_id==neuron_id) {
	multimeter_id = header_multi[thread_num].nodes[i].owner_id;
	numberOfValues = header_multi[thread_num].nodes[i].numberOfValues;
      }
    }
  if (!buffer_multi[thread_num].isEnoughFreeSpace(4*sizeof(int)+numberOfValues*sizeof(double)))
  {
    sion_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
    buffer_multi[thread_num].clear();
  }
  buffer_multi[thread_num] << multimeter_id << neuron_id << timestamp << numberOfValues;
  for (int i=0; i<numberOfValues; i++)
  {
    buffer_multi[thread_num] << v[i];
  }
}


void Sionlib_logger::record_spike(int neuron_id, double t)
{
  const int thread_num = omp_get_thread_num();
    switch (loggerType) {
      case nestio::Standard:
	srecord_spike(neuron_id, t);
	break;
      case nestio::Buffered:
	brecord_spike(neuron_id, t);
	break;
      case nestio::Collective:
	crecord_spike(neuron_id, t);
	break;
    }
    header_spike[thread_num].numberOfWrittenData++;
}

void Sionlib_logger::record_multi(int neuron_id, int timestamp, double* v)
{
  const int thread_num = omp_get_thread_num();
    switch (loggerType) {
      case nestio::Standard:
	srecord_multi(neuron_id, timestamp, v);
	break;
      case nestio::Buffered:
	brecord_multi(neuron_id, timestamp, v);
	break;
      case nestio::Collective:
	crecord_multi(neuron_id, timestamp, v);
	break;
    }
    header_multi[thread_num].numberOfWrittenData++;
}

void Sionlib_logger::srecord_spike(int neuron_id, double t)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif

  sion_ensure_free_space(spike_sid[thread_num], 2*sizeof(int)+sizeof(double));
  int spikedetector_id = -1;
  for (int i=0; i<header_spike[thread_num].NodesCount;i++) {
    if (header_spike[thread_num].nodes[i].owner_id==neuron_id) {
      spikedetector_id = header_spike[thread_num].nodes[i].owner_id;
    }
  }
  sion_fwrite(&spikedetector_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&t, sizeof(double), 1, spike_sid[thread_num]);
}

void Sionlib_logger::srecord_multi(int neuron_id, int timestamp, double* v)
{
  #ifdef _OPENMP
  const int thread_num = omp_get_thread_num();
  #else
  const int thread_num = 0;
  #endif
  int multimeter_id=-1;
  int numberOfValues=0;
  for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
    if (header_multi[thread_num].nodes[i].owner_id==neuron_id) {
      multimeter_id = header_multi[thread_num].nodes[i].owner_id;
      numberOfValues = header_multi[thread_num].nodes[i].numberOfValues;
    }
  }
  
  sion_ensure_free_space(multi_sid[thread_num], 4*sizeof(int)+(numberOfValues)*sizeof(double));
  sion_fwrite(&multimeter_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&timestamp, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&numberOfValues, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(v, sizeof(double), numberOfValues, multi_sid[thread_num]);
  
  header_multi[thread_num].numberOfWrittenData++;
}

/*
 * Collecting header informations for SpikeDetectors
 */
void Sionlib_logger::signup_spike(SpikeDetector* spike, int neuron_id, int buf)
{
  const int thread_num = omp_get_thread_num();
#pragma omp critical
  {
    header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].owner_id=neuron_id;
    header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].id=spike->spikedetector_id;
    header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].interval=0;
    header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].type=2;
    header_spike[thread_num].NodesCount++;
  }
}

/*
 * Collecting header informations for Multimeters
 */
void Sionlib_logger::signup_multi(Multimeter* multi, int neuron_id, int buf)
{
  const int thread_num = omp_get_thread_num();
#pragma omp critical
  {
    header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].id=multi->multimeter_id;
    header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].owner_id=neuron_id;
    header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].interval=multi->samlpingInterval;
    header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].numberOfValues=multi->numberOfValues;   //Work around
    for (int i=0; i<multi->numberOfValues; i++) {
      memcpy(header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].valueNames[i],multi->valueNames.at(i).c_str(), min(20,(int)multi->valueNames.at(i).size()));
    }
    header_multi[thread_num].NodesCount++;
  }
}

/*
 *  Write header to sion files: Spikedetector and Multimeter file
 */
void Sionlib_logger::createDatasets()
{
  
  const int thread_num = omp_get_thread_num();
  
  std::cout << "createDatasets" << std::endl;
  int numberOfRecords=-1;
  int startOfBody;
  //
  // create Spike createDatasets
  //
  sion_fwrite(&header_spike[thread_num].NodesCount, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&simSettings.T, sizeof(double), 1, spike_sid[thread_num]);
  sion_fwrite(&simSettings.Tresolution, sizeof(double), 1, spike_sid[thread_num]); // should be Tresolution
  sion_fwrite(&numberOfRecords, sizeof(int), 1, spike_sid[thread_num]);
  
  // !!!!! caution
  startOfBody = 3*sizeof(int)+2*sizeof(double)+header_spike[thread_num].NodesCount*(2*sizeof(int));
  sion_fwrite(&startOfBody, sizeof(int), 1, spike_sid[thread_num]);
  
  /*for (int i=0; i<header_spike.NodesCount;i++) {
    sion_fwrite(&header_spike.nodes[i].id, sizeof(int), 1, spike_sid);
    sion_fwrite(&header_spike.nodes[i].owner_id, sizeof(int), 1, spike_sid);
  }*/
  
  //
  // create Multidatasets
  //
  
  sion_fwrite(&header_multi[thread_num].NodesCount, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&simSettings.T, sizeof(double), 1, multi_sid[thread_num]);
  sion_fwrite(&simSettings.Tresolution, sizeof(double), 1, multi_sid[thread_num]); // should be Tresolution
  sion_fwrite(&numberOfRecords, sizeof(int), 1, multi_sid[thread_num]);
  
  // !!!!! caution startOfBody has to be set correctly
  startOfBody = 3*sizeof(int)+2*sizeof(double)+header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double));
  for (int i=0; i<header_multi[thread_num].NodesCount;i++)
    startOfBody += header_multi[thread_num].NodesCount*header_multi[thread_num].nodes[i].numberOfValues*sizeof(char)*20;
  
  sion_fwrite(&startOfBody, sizeof(int), 1, multi_sid[thread_num]);
  
  for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
    sion_fwrite(&header_multi[thread_num].nodes[i].id, sizeof(int), 1, multi_sid[thread_num]);
    //sion_fwrite(&header_multi.nodes[i].owner_id, sizeof(int), 1, multi_sid);
    sion_fwrite(&header_multi[thread_num].nodes[i].interval, sizeof(double), 1, multi_sid[thread_num]);
    sion_fwrite(&header_multi[thread_num].nodes[i].numberOfValues, sizeof(int), 1, multi_sid[thread_num]);
    for (int j=0; j<header_multi[thread_num].nodes[i].numberOfValues; j++) {
      sion_fwrite(header_multi[thread_num].nodes[i].valueNames[j], sizeof(char), 20, multi_sid[thread_num]);
    }
      
  }
}

/*
 * Init sion file
 */
Sionlib_logger::Sionlib_logger(std::string filename, int ibuf_size, nestio::SionLoggerType loggerType, nestio::SimSettings &simSettings)
:simSettings(simSettings), loggerType(loggerType)
{
	//std::cout << "create logger in thread "<< omp_get_thread_num() << std::endl;
	int num_threads=omp_get_max_threads();
    
	header_multi.resize(num_threads);
	header_spike.resize(num_threads);
	
	for (int i=0; i<num_threads; i++) {
	  header_multi[i].NodesCount=0;
	  header_spike[i].NodesCount=0;
	  
	  header_multi[i].numberOfWrittenData=0;
	  header_spike[i].numberOfWrittenData=0;
	}

	spike_sid.reserve(num_threads);
	multi_sid.reserve(num_threads);
	
	#pragma omp parallel
	{
	
	  int thread_num = omp_get_thread_num();
	  
	  
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
	  numFiles = 1; //internal sion can write to more than one file -> filesystem/performance
	  fsblksize = -1;
	  
	  std::stringstream spike_ss, multi_ss;
	  spike_ss << "spikes_" << filename;
	  multi_ss << "multi_" << filename;
	  
	  char spike_fname[256],multi_fname[256];
	  strcpy(spike_fname, spike_ss.str().c_str());
	  strcpy(multi_fname, multi_ss.str().c_str());
	  
	  /* create a new file */
	  spike_sid[thread_num] = sion_paropen_ompi(spike_fname, "bw", &numFiles,
	  MPI_COMM_WORLD, &lComm,
	  &buf_size, &fsblksize,
	  &own_id,
	  NULL, &newfname);
	  
	  /* create a new file */
	  multi_sid[thread_num] = sion_paropen_ompi(multi_fname, "bw", &numFiles,
	  MPI_COMM_WORLD, &lComm,
	  &buf_size, &fsblksize,
	  &own_id,
	  NULL, &newfname);
	}
	
	if (loggerType == nestio::Buffered || loggerType == nestio::Collective) {
	  buffer_multi = new SionBuffer[num_threads];
	  buffer_spike = new SionBuffer[num_threads];
	  for (int i=0; i<num_threads;i++) {
	    buffer_multi[i].extend(buf_size);
	    buffer_spike[i].extend(buf_size);
	  }
	}
}


/*
 * close sion files
 */
Sionlib_logger::~Sionlib_logger()
{
  #pragma omp parallel
  {
    const int thread_num = omp_get_thread_num();
    
    //write values from buffer to file
    sion_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
    buffer_spike[thread_num].clear();
    
    sion_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
    buffer_multi[thread_num].clear();
    
    //write TAIL
    sion_fwrite(&header_spike[thread_num].numberOfWrittenData,sizeof(int),1,spike_sid[thread_num]);
    sion_fwrite(&header_multi[thread_num].numberOfWrittenData,sizeof(int),1,multi_sid[thread_num]);
    
    sion_parclose_ompi(multi_sid[thread_num]);
    sion_parclose_ompi(spike_sid[thread_num]);
  }
}


/*
 * called during the sync step
 * could be used to call sion_ensure_free_space
 */
void Sionlib_logger::updateDatasetSizes(const double& t)
{
  if (loggerType == nestio::Collective) {
    #ifdef _OPENMP
    const int thread_num = omp_get_thread_num();
    #else
    const int thread_num = 0;
    #endif
    
    sion_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
    buffer_spike[thread_num].clear();
    
    sion_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
    buffer_multi[thread_num].clear();
    
  }
}

std::ostream& operator << (std::ostream &o, const Sionlib_logger &l)
{
  o << "Sionlib_logger";
  return o;
}
