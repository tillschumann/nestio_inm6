#include "sionlib_logger.h"
#include "iostream"
#include "mpi.h"
#include "sstream"
#include <omp.h>


template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

void Sionlib_logger::crecord_spike(int spikedetector_id, int neuron_id, int timestamp)
{
  const int thread_num = omp_get_thread_num();
 
  buffer_spike[thread_num].getEnoughFreeSpace(3*sizeof(int));
  buffer_spike[thread_num] << spikedetector_id << neuron_id << timestamp;
}

void Sionlib_logger::brecord_spike(int spikedetector_id, int neuron_id, int timestamp)
{
  const int thread_num = omp_get_thread_num();
  
  if (!buffer_spike[thread_num].isEnoughFreeSpace(3*sizeof(int)))
  {
    sion_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
    buffer_spike[thread_num].clear();
  }
  buffer_spike[thread_num] << spikedetector_id << neuron_id << timestamp;
}

void Sionlib_logger::crecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  const int thread_num = omp_get_thread_num();
  
  int numberOfValues=data.size();
  buffer_multi[thread_num].getEnoughFreeSpace(4*sizeof(int)+numberOfValues*sizeof(double));
  buffer_multi[thread_num] << multimeter_id << neuron_id << timestamp << numberOfValues;
  for (int i=0; i<numberOfValues; i++)
  {
    buffer_multi[thread_num] << data[i];
  }
}

void Sionlib_logger::brecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  const int thread_num = omp_get_thread_num();
  
  int numberOfValues=data.size();
  if (!buffer_multi[thread_num].isEnoughFreeSpace(4*sizeof(int)+numberOfValues*sizeof(double)))
  {
    sion_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
    buffer_multi[thread_num].clear();
  }
  buffer_multi[thread_num] << multimeter_id << neuron_id << timestamp << numberOfValues;
  for (int i=0; i<numberOfValues; i++)
  {
    buffer_multi[thread_num] << data[i];
  }
}


void Sionlib_logger::record_spike(int spikedetector_id, int neuron_id, int timestamp)
{
  const int thread_num = omp_get_thread_num();
    switch (P_.loggerType_) {
      case nestio::Standard:
	srecord_spike(spikedetector_id, neuron_id, timestamp);
	break;
      case nestio::Buffered:
	brecord_spike(spikedetector_id, neuron_id, timestamp);
	break;
      case nestio::Collective:
	crecord_spike(spikedetector_id, neuron_id, timestamp);
	break;
    }
    header_spike[thread_num].numberOfWrittenData++;
}

void Sionlib_logger::record_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  const int thread_num = omp_get_thread_num();
    switch (P_.loggerType_) {
      case nestio::Standard:
	srecord_multi(multimeter_id,neuron_id, timestamp, data);
	break;
      case nestio::Buffered:
	brecord_multi(multimeter_id,neuron_id, timestamp, data);
	break;
      case nestio::Collective:
	crecord_multi(multimeter_id,neuron_id, timestamp, data);
	break;
    }
    header_multi[thread_num].numberOfWrittenData++;
  
}

void Sionlib_logger::srecord_spike(int spikedetector_id, int neuron_id, int timestamp)
{
  const int thread_num = omp_get_thread_num();

  sion_ensure_free_space(spike_sid[thread_num], 3*sizeof(int));
  sion_fwrite(&spikedetector_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&timestamp, sizeof(int), 1, spike_sid[thread_num]);
}

void Sionlib_logger::srecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  const int thread_num = omp_get_thread_num();

  int numberOfValues=data.size();
  
  sion_ensure_free_space(multi_sid[thread_num], 4*sizeof(int)+(numberOfValues)*sizeof(double));
  sion_fwrite(&multimeter_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&timestamp, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&numberOfValues, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&data[0], sizeof(double), numberOfValues, multi_sid[thread_num]);
}

/*
 * Collecting header informations for SpikeDetectors
 */
void Sionlib_logger::signup_spike(int id, int neuron_id, int expectedSpikeCount)
{
  //std::cout << "Sionlib_logger::signup_spike" << std::endl;
  const int thread_num = omp_get_thread_num();
#pragma omp critical
  {
    SionFileHeaderNode node;
    node.neuron_id=neuron_id;
    node.id=id;
    node.interval=0;
    
    header_spike[thread_num].nodes.push_back(node);
    
    
    //header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].neuron_id=neuron_id;
    //header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].id=spike->spikedetector_id;
    //header_spike[thread_num].nodes[header_spike[thread_num].NodesCount].interval=0;
    header_spike[thread_num].NodesCount++;
  }
}

/*
 * Collecting header informations for Multimeters
 */
void Sionlib_logger::signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames, double simulationTime)
{
  const int thread_num = omp_get_thread_num();
#pragma omp critical
  {
    SionFileHeaderNode node;
    node.id=id;
    node.neuron_id=neuron_id;
    node.interval=sampling_interval;
    node.numberOfValues=valueNames.size();
    for (int i=0; i<valueNames.size(); i++) {
      memcpy(node.valueNames[i],valueNames.at(i).c_str(), min(20,(int)valueNames.at(i).size()));
    }
    
    header_multi[thread_num].nodes.push_back(node);
    
    
    //header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].id=multi->multimeter_id;
    //header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].neuron_id=neuron_id;
    //header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].interval=multi->samlpingInterval;
    //header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].numberOfValues=multi->numberOfValues;   //Work around
    //for (int i=0; i<multi->numberOfValues; i++) {
    //  memcpy(header_multi[thread_num].nodes[header_multi[thread_num].NodesCount].valueNames[i],multi->valueNames.at(i).c_str(), min(20,(int)multi->valueNames.at(i).size()));
    //}
    header_multi[thread_num].NodesCount++;
  }
}

/*
 *  Write header to sion files: Spikedetector and Multimeter file
 */
void Sionlib_logger::initialize(const double T)
{
  #pragma omp parallel
  {
  
    int thread_num = omp_get_thread_num();
    
    
    /* SION parameters */
    int numFiles, rank, num_procs;
    MPI_Comm lComm;
    //sion_int64 left, bwrote;
    sion_int32 fsblksize;
    char *newfname=NULL;
    /* MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    /* open parameters */
    numFiles = 1; //internal sion can write to more than one file -> filesystem/performance
    fsblksize = -1;
    
    
    std::string sfn;
    std::string mfn;
    
    char spike_fname[256],multi_fname[256];
    strcpy(spike_fname, sfn.c_str());
    strcpy(multi_fname, mfn.c_str());
    
    /* create a new file */
    spike_sid[thread_num] = sion_paropen_ompi(spike_fname, "bw", &numFiles,
    MPI_COMM_WORLD, &lComm,
    &P_.sion_buffer_size_, &fsblksize,
    &rank,
    NULL, &newfname);
    
    /* create a new file */
    multi_sid[thread_num] = sion_paropen_ompi(multi_fname, "bw", &numFiles,
    MPI_COMM_WORLD, &lComm,
    &P_.sion_buffer_size_, &fsblksize,
    &rank,
    NULL, &newfname);

    #ifdef _DEBUG_MODE
    std::cout << "createDatasets" << std::endl;
    #endif
    int numberOfRecords=-1;
    int startOfBody;
    //
    // create Spike createDatasets
    //
    sion_fwrite(&header_spike[thread_num].NodesCount, sizeof(int), 1, spike_sid[thread_num]);
    sion_fwrite(&P_.T_, sizeof(double), 1, spike_sid[thread_num]);
    sion_fwrite(&P_.Tresolution_, sizeof(double), 1, spike_sid[thread_num]); // should be Tresolution //TODO
    sion_fwrite(&numberOfRecords, sizeof(int), 1, spike_sid[thread_num]);

    // !!!!! caution
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_spike[thread_num].NodesCount*(2*sizeof(int));
    sion_fwrite(&startOfBody, sizeof(int), 1, spike_sid[thread_num]);

    /*for (int i=0; i<header_spike.NodesCount;i++) {
      sion_fwrite(&header_spike.nodes[i].id, sizeof(int), 1, spike_sid);
      sion_fwrite(&header_spike.nodes[i].neuron_id, sizeof(int), 1, spike_sid);
    }*/

    //
    // create Multidatasets
    //

    sion_fwrite(&header_multi[thread_num].NodesCount, sizeof(int), 1, multi_sid[thread_num]);
    sion_fwrite(&P_.T_, sizeof(double), 1, multi_sid[thread_num]);
    sion_fwrite(&P_.Tresolution_, sizeof(double), 1, multi_sid[thread_num]); // should be Tresolution
    sion_fwrite(&numberOfRecords, sizeof(int), 1, multi_sid[thread_num]);

    // !!!!! caution startOfBody has to be set correctly
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double));
    for (int i=0; i<header_multi[thread_num].NodesCount;i++)
      startOfBody += header_multi[thread_num].NodesCount*header_multi[thread_num].nodes[i].numberOfValues*sizeof(char)*20;

    sion_fwrite(&startOfBody, sizeof(int), 1, multi_sid[thread_num]);

    for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
      sion_fwrite(&header_multi[thread_num].nodes[i].id, sizeof(int), 1, multi_sid[thread_num]);
      sion_fwrite(&header_multi[thread_num].nodes[i].neuron_id, sizeof(int), 1, multi_sid[thread_num]);
      sion_fwrite(&header_multi[thread_num].nodes[i].interval, sizeof(double), 1, multi_sid[thread_num]);
      sion_fwrite(&header_multi[thread_num].nodes[i].numberOfValues, sizeof(int), 1, multi_sid[thread_num]);
      for (int j=0; j<header_multi[thread_num].nodes[i].numberOfValues; j++) {
	sion_fwrite(header_multi[thread_num].nodes[i].valueNames[j], sizeof(char), 20, multi_sid[thread_num]);
      }
	
    }
  }
}

/*
 * Init sion file
 */
Sionlib_logger::Sionlib_logger(const std::string& path, const std::string& file_extension,int logger_buf_size, sion_int64 sion_buf_size, nestio::LoggerType loggerType)
:P_(path, file_extension, loggerType, logger_buf_size, sion_buf_size)
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
    
    if (P_.loggerType_ == nestio::Buffered || P_.loggerType_ == nestio::Collective) {
    buffer_multi = new SionBuffer[num_threads];
    buffer_spike = new SionBuffer[num_threads];
    for (int i=0; i<num_threads;i++) {
      buffer_multi[i].extend(P_.logger_buffer_size_);
      buffer_spike[i].extend(P_.logger_buffer_size_);
    }
  }
}


/*
 * close sion files
 */
Sionlib_logger::~Sionlib_logger()
{
  if (P_.loggerType_ == nestio::Buffered || P_.loggerType_ == nestio::Collective) {
    delete buffer_multi;
    delete buffer_spike;
  }
}


void Sionlib_logger::finalize()
{
  #pragma omp parallel
  {
    const int thread_num = omp_get_thread_num();
    
    if (P_.loggerType_ == nestio::Collective) {
      sion_coll_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
      buffer_spike[thread_num].clear();
      sion_coll_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
      buffer_multi[thread_num].clear();
    }
    else if (P_.loggerType_ == nestio::Buffered) {
      //write values from buffer to file 
      if (buffer_spike[thread_num].getSize()>0) {   
	sion_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
	buffer_spike[thread_num].clear();
      }

      if (buffer_multi[thread_num].getSize()>0) {
	sion_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
	buffer_multi[thread_num].clear();
      }
    }
    
    //write TAIL
    #ifdef _DEBUG_MODE
    std::cout << "spike numberOfWrittenData=" << header_spike[thread_num].numberOfWrittenData << std::endl;
    std::cout << "multi numberOfWrittenData=" << header_multi[thread_num].numberOfWrittenData << std::endl;
    #endif
    
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
void Sionlib_logger::syncronize(const double& t)
{
  if (P_.loggerType_ == nestio::Collective) {
    const int thread_num = omp_get_thread_num();

    //if (buffer_spike[thread_num].getSize()>0) {
    sion_coll_fwrite(buffer_spike[thread_num].read(), buffer_spike[thread_num].getSize(), 1, spike_sid[thread_num]);
    buffer_spike[thread_num].clear();
    //}
    //if (buffer_multi[thread_num].getSize()>0) {
    sion_coll_fwrite(buffer_multi[thread_num].read(), buffer_multi[thread_num].getSize(), 1, multi_sid[thread_num]);
    buffer_multi[thread_num].clear();
    //}
    
  }
}

std::ostream& operator << (std::ostream &o, const Sionlib_logger &l)
{
  o << "Sionlib_logger";
  return o;
}

Sionlib_logger::Parameters_::Parameters_(const std::string& path, const std::string& file_extension, nestio::LoggerType loggerType, int logger_buffer_size, sion_int64 sion_buffer_size)
: loggerType_(loggerType), path_(path), file_extension_(file_extension), logger_buffer_size_(logger_buffer_size), sion_buffer_size_(sion_buffer_size)
{}
