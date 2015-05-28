#include "sionlib_logger.h"

#ifdef HAVE_SIONLIB
#include "iostream"
#include "sstream"
#include <omp.h>

#ifndef NESTIOPROXY
#include "dictutils.h"
#include "network.h"
#include "node.h"
#endif


#ifdef HAVE_MPI
#include <mpi.h>
#endif /* #ifdef HAVE_MPI */


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

  //sion_ensure_free_space(spike_sid[thread_num], 3*sizeof(int));
  sion_fwrite(&spikedetector_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, spike_sid[thread_num]);
  sion_fwrite(&timestamp, sizeof(int), 1, spike_sid[thread_num]);
}

void Sionlib_logger::srecord_multi(int multimeter_id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  const int thread_num = omp_get_thread_num();

  int numberOfValues=data.size();
  
  //sion_ensure_free_space(multi_sid[thread_num], 4*sizeof(int)+(numberOfValues)*sizeof(double));
  sion_fwrite(&multimeter_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&neuron_id, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&timestamp, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&numberOfValues, sizeof(int), 1, multi_sid[thread_num]);
  sion_fwrite(&data[0], sizeof(double), numberOfValues, multi_sid[thread_num]);
}

/*
 * Collecting header informations for SpikeDetectors
 */
void Sionlib_logger::signup_spike(int id, int neuron_id)
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
    
   
    header_spike[thread_num].NodesCount++;
  }
}

/*
 * Collecting header informations for Multimeters
 */
void Sionlib_logger::signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames)
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
      #ifndef NESTIOPROXY
      std::string vN = valueNames.at(i).toString(); //TODO
      #else
      std::string vN = valueNames.at(i);		//TODO
      #endif
      memcpy(node.valueNames[i],vN.c_str(), min(20,(int)vN.size())); //TODO
    }
    
    header_multi[thread_num].nodes.push_back(node);
    
    header_multi[thread_num].NodesCount++;
  }
}

/*
 *  write header to sion files: Spikedetector and Multimeter file
 */
void Sionlib_logger::writeHeaders2File(const int& thread_num)
{
  #ifdef _DEBUG_MODE
  std::cout << "writeHead2File" << std::endl;
  #endif
  int numberOfRecords=-1;
  int startOfBody;
  //
  // create Spike createDatasets
  //	write headers for spikedetectors
  #ifdef _SIONLIB_COLL
  if (P_.loggerType_ == nestio::Collective) {
    sion_coll_fwrite(&header_spike[thread_num].NodesCount, sizeof(int), 1, spike_sid[thread_num]);
    sion_coll_fwrite(&P_.T_, sizeof(double), 1, spike_sid[thread_num]);
    sion_coll_fwrite(&P_.Tresolution_, sizeof(double), 1, spike_sid[thread_num]); 
    sion_coll_fwrite(&numberOfRecords, sizeof(int), 1, spike_sid[thread_num]);

    // !!!!! caution
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_spike[thread_num].NodesCount*(2*sizeof(int));
    sion_coll_fwrite(&startOfBody, sizeof(int), 1, spike_sid[thread_num]);
  }
  else
  #endif
  {
    sion_fwrite(&header_spike[thread_num].NodesCount, sizeof(int), 1, spike_sid[thread_num]);
    sion_fwrite(&P_.T_, sizeof(double), 1, spike_sid[thread_num]);
    sion_fwrite(&P_.Tresolution_, sizeof(double), 1, spike_sid[thread_num]);
    sion_fwrite(&numberOfRecords, sizeof(int), 1, spike_sid[thread_num]);

    // !!!!! caution
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_spike[thread_num].NodesCount*(2*sizeof(int));
    sion_fwrite(&startOfBody, sizeof(int), 1, spike_sid[thread_num]);
  }

  //
  // create Multidatasets
  // write headers for multimeters
  
  #ifdef _SIONLIB_COLL
  if (P_.loggerType_ == nestio::Collective) {
    sion_coll_fwrite(&header_multi[thread_num].NodesCount, sizeof(int), 1, multi_sid[thread_num]);
    sion_coll_fwrite(&P_.T_, sizeof(double), 1, multi_sid[thread_num]);
    sion_coll_fwrite(&P_.Tresolution_, sizeof(double), 1, multi_sid[thread_num]); // should be Tresolution
    sion_coll_fwrite(&numberOfRecords, sizeof(int), 1, multi_sid[thread_num]);

    // !!!!! caution startOfBody has to be set correctly
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double));
    for (int i=0; i<header_multi[thread_num].NodesCount;i++)
      startOfBody += header_multi[thread_num].NodesCount*header_multi[thread_num].nodes[i].numberOfValues*sizeof(char)*20;

    sion_coll_fwrite(&startOfBody, sizeof(int), 1, multi_sid[thread_num]);

    SionBuffer buffer(header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double)));

    for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
      SionFileHeaderNode& node = header_multi[thread_num].nodes[i];
      buffer.getEnoughFreeSpace(3*sizeof(int)+sizeof(double)+sizeof(char)*node.numberOfValues);
      buffer << node.id;
      buffer << node.neuron_id;
      buffer << node.interval;
      buffer << node.numberOfValues;
      for (int j=0; j<node.numberOfValues; j++)
	buffer << node.valueNames[j];
    }
    if (buffer.getSize()>0)
      sion_coll_fwrite(buffer.read(), buffer.getSize(), 1, multi_sid[thread_num]);
    else
      sion_coll_fwrite(buffer.read(), 1, 0, multi_sid[thread_num]);
  }
  else
  #endif
  {
    sion_fwrite(&header_multi[thread_num].NodesCount, sizeof(int), 1, multi_sid[thread_num]);
    sion_fwrite(&P_.T_, sizeof(double), 1, multi_sid[thread_num]);
    sion_fwrite(&P_.Tresolution_, sizeof(double), 1, multi_sid[thread_num]); // should be Tresolution
    sion_fwrite(&numberOfRecords, sizeof(int), 1, multi_sid[thread_num]);

    // !!!!! caution startOfBody has to be set correctly
    startOfBody = 3*sizeof(int)+2*sizeof(double)+header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double));
    for (int i=0; i<header_multi[thread_num].NodesCount;i++)
      startOfBody += header_multi[thread_num].NodesCount*header_multi[thread_num].nodes[i].numberOfValues*sizeof(char)*20;

    sion_fwrite(&startOfBody, sizeof(int), 1, multi_sid[thread_num]);

    SionBuffer buffer(header_multi[thread_num].NodesCount*(3*sizeof(int)+sizeof(double)));
    //buffer << (int)1;
    for (int i=0; i<header_multi[thread_num].NodesCount;i++) {
      SionFileHeaderNode& node = header_multi[thread_num].nodes[i];
      buffer.getEnoughFreeSpace(3*sizeof(int)+sizeof(double)+20*sizeof(char)*node.numberOfValues);
      
      buffer << node.id;
      buffer << node.neuron_id;
      buffer << node.interval;
      buffer << node.numberOfValues;
      for (int j=0; j<node.numberOfValues; j++)
	buffer << node.valueNames[j];
    }
    //std::cout << "tricky size=" << buffer.getSize() << std::endl;
    if (buffer.getSize()>0)
      sion_fwrite(buffer.read(), buffer.getSize(), 1, multi_sid[thread_num]);
    else
      sion_fwrite(buffer.read(), 1, 0, multi_sid[thread_num]);
  }
}

/*
 *  Open Sion file and
 *  write header to sion files: Spikedetector and Multimeter file
 */
void Sionlib_logger::initialize(const double T)
{  
  //set missing parameter
  P_.T_ = T;
  
  
  #pragma omp parallel
  {
  
    int thread_num = omp_get_thread_num();
    
    
    //generate file names
    char spike_fname[256],multi_fname[256];
    strcpy(spike_fname, build_filename_("spikedetector").c_str());
    strcpy(multi_fname, build_filename_("multimeter").c_str());
    
    #ifdef _DEBUG_MODE
    std::cout << "Sionlib_logger::initialize" << std::endl;
    std::cout << "filepath spike: " << spike_fname << std::endl;
    std::cout << "filepath multi: " << multi_fname << std::endl;
    #endif
    
    
    
    
    /* SION parameters */
        
    //sion_int64 left, bwrote;
    sion_int32 fsblksize;
    char *newfname=NULL;
    
  
    
    /* open parameters */
    int numFiles = 1; //internal sion can write to more than one file -> filesystem/performance
    fsblksize = -1;
    
    /* MPI */
    #ifdef HAVE_MPI
    int rank, num_procs;
    
    sion_int64 sion_buffer_size_spike = P_.sion_buffer_size_;
    sion_int64 sion_buffer_size_multi = P_.sion_buffer_size_;
    
    MPI_Comm lComm;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    /* create a new file */
    spike_sid[thread_num] = sion_paropen_ompi(spike_fname, "bw", &numFiles,
    MPI_COMM_WORLD, &lComm,
    &sion_buffer_size_spike, &fsblksize,
    &rank,
    NULL, &newfname);
    
    #ifdef _DEBUG_MODE
    std::cout << "Sionlib_logger::initialize created spike_sid" << std::endl;
    #endif
    
    /* create a new file */
    multi_sid[thread_num] = sion_paropen_ompi(multi_fname, "bw", &numFiles,
    MPI_COMM_WORLD, &lComm,
    &sion_buffer_size_multi, &fsblksize,
    &rank,
    NULL, &newfname);
    
    
    #ifdef _DEBUG_MODE
    std::cout << "Sionlib_logger::initialize created multi_sid" << std::endl;
    #endif
    
    #else
    int ntasks=1;
    int rank = 0;
    int* ptr_rank = &rank;
    
    sion_int64* ptr_chunk = &P_.sion_buffer_size_;
    
    spike_sid[thread_num] =  sion_open(spike_fname, "bw", &ntasks, &numFiles,
    &ptr_chunk, &fsblksize,
    &ptr_rank,
    NULL);
    
    multi_sid[thread_num] = sion_open(multi_fname, "bw", &ntasks,  &numFiles,
    &ptr_chunk, &fsblksize,
    &ptr_rank,
    NULL);
    
    #endif

    //write headers to files
    writeHeaders2File(thread_num);
  }
}

/*
 * Init sion file
 */
Sionlib_logger::Sionlib_logger(const std::string& path, const std::string& file_extension,int logger_buf_size, sion_int64 sion_buf_size, nestio::Logger_type loggerType)
:P_(path, file_extension, loggerType, logger_buf_size, sion_buf_size)
{  
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

Sionlib_logger::Sionlib_logger(): P_(".",".log", nestio::Standard, 100, 100) {
  //std::cout << "create logger in thread "<< omp_get_thread_num() << std::endl;
  
}


/*
 * close sion files
 */
Sionlib_logger::~Sionlib_logger()
{
  if (P_.loggerType_ == nestio::Buffered || P_.loggerType_ == nestio::Collective) {
    delete [] buffer_multi;
    delete [] buffer_spike;
  }
}

/*
 * write buffers to file and close files
 * 
 */
void Sionlib_logger::finalize()
{
  #pragma omp parallel
  {
    const int thread_num = omp_get_thread_num();
    #ifdef _SIONLIB_COLL
    if (P_.loggerType_ == nestio::Collective) {
      sion_coll_fwrite(buffer_spike[thread_num].read(), 1, buffer_spike[thread_num].getSize(), spike_sid[thread_num]);
      buffer_spike[thread_num].clear();
      sion_coll_fwrite(buffer_multi[thread_num].read(), 1, buffer_multi[thread_num].getSize(), multi_sid[thread_num]);
      buffer_multi[thread_num].clear();
    }
    else
    #endif  
    if (P_.loggerType_ == nestio::Buffered) {
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
    #ifdef _SIONLIB_COLL
    if (P_.loggerType_ == nestio::Collective) {
      sion_coll_fwrite(&header_spike[thread_num].numberOfWrittenData,sizeof(int),1,spike_sid[thread_num]);
      sion_coll_fwrite(&header_multi[thread_num].numberOfWrittenData,sizeof(int),1,multi_sid[thread_num]);
    }
    else
    #endif
    {
      sion_fwrite(&header_spike[thread_num].numberOfWrittenData,sizeof(int),1,spike_sid[thread_num]);
      sion_fwrite(&header_multi[thread_num].numberOfWrittenData,sizeof(int),1,multi_sid[thread_num]);
    }
    #ifdef HAVE_MPI
    sion_parclose_ompi(multi_sid[thread_num]);
    sion_parclose_ompi(spike_sid[thread_num]);
    #else
    sion_close(multi_sid[thread_num]);
    sion_close(spike_sid[thread_num]);
    #endif
  }
}


/*
 * called during the sync step
 * could be used to call sion_ensure_free_space
 * for collective: writes buffer to disk 
 * 
 */
void Sionlib_logger::synchronize(const double& t)
{
  #ifdef _DEBUG_MODE
  std::cout << "synchronize t=" << t << std::endl;
  #endif
  
  #ifdef _SIONLIB_COLL
  if (P_.loggerType_ == nestio::Collective) {
    const int thread_num = omp_get_thread_num();

    if (buffer_spike[thread_num].getSize()>0) {
      sion_coll_fwrite(buffer_spike[thread_num].read(), 1, buffer_spike[thread_num].getSize(), spike_sid[thread_num]);
      buffer_spike[thread_num].clear();
    }
    if (buffer_multi[thread_num].getSize()>0) {
      sion_coll_fwrite(buffer_multi[thread_num].read(), 1, buffer_multi[thread_num].getSize(), multi_sid[thread_num]);
      buffer_multi[thread_num].clear();
    }
  }
  #endif
}

std::ostream& operator << (std::ostream &o, const Sionlib_logger &l)
{
  o << "Sionlib_logger";
  return o;
}

#ifndef NESTIOPROXY
void Sionlib_logger::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;    // temporary copy in case of errors
  ptmp.set(d);   // throws if BadProperty
  
  P_ = ptmp;
}
#endif

Sionlib_logger::Parameters_::Parameters_(const std::string& path, const std::string& file_extension, nestio::Logger_type loggerType, int logger_buffer_size, sion_int64 sion_buffer_size)
: loggerType_(loggerType), path_(path), file_extension_(file_extension), logger_buffer_size_(logger_buffer_size), sion_buffer_size_(sion_buffer_size)
{}


 #ifndef NESTIOPROXY	
void Sionlib_logger::Parameters_::set(const DictionaryDatum& d)
{
  std::string loggerType_str;
  updateValue<std::string>(d, "buf", loggerType_str);
  if (loggerType_str == "Bufferd")
      loggerType_ = nest::Buffered;
  #ifdef _SIONLIB_COLL
  else if (loggerType_str == "Collective")
      loggerType_ = nestio::Collective;
  #endif
  else
      loggerType_ = nestio::Standard;
  
  //get resolution from network
  P_.Tresolution_ = nest::Node::network()->get_min_delay();
  
  updateValue<bool>(d, "overwrite_files", overwrite_files_);
  updateValue<std::string>(d, "path", path_);
  updateValue<std::string>(d, "file_extension", file_extension_);
  long sion_buffer_size_long;
  updateValue<long>(d, "sion_buffer_size", sion_buffer_size_long);
  sion_buffer_size_ = sion_buffer_size_long;
  #ifdef _SIONLIB_COLL
  if (loggerType_ == nestio::Collective || loggerType_ == nestio::Buffered) {
  #else
  if (loggerType_ == nestio::Buffered) {
  #endif
    long logger_buffer_size_long;
    updateValue<long>(d, "logger_buffer_size", logger_buffer_size_long);
    logger_buffer_size_ = logger_buffer_size_long;
    
  }
  
  //TODO error handling
}
#endif

const std::string Sionlib_logger::build_filename_(std::string prefix) const
{
  // number of digits in number of virtual processes
#ifndef NESTIOPROXY
  
  std::ostringstream basename;
  const std::string& path = Node::network()->get_data_path();
  if (!P_.path_.empty())
    basename << P_.path_ << '/';
  else if ( !path.empty() )
    basename << path << '/';
  basename << prefix;

  return basename.str() + '.' + P_.file_extension_;
#else
  
  std::ostringstream file_path;
  file_path << P_.path_ << '/';
  file_path << prefix << "." << P_.file_extension_;
  return file_path.str();
#endif
}

#endif /* #ifdef HAVE_SIONLIB */
