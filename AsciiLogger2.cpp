#include "AsciiLogger2.h"
#include "iostream"
#include "sstream"
#include <fstream> 
#include <iomanip>
#ifndef NESTIOPROXY
#include "network.h"
#endif
//#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

void nest::AsciiLogger2::record_spike(int id, int neuron_id, int timestamp)
{

  std::ofstream& fs_ = *(spikedetectors[id].fs_);
  #pragma omp critical
  {
    fs_ << id;
    fs_ << " " << neuron_id;
    fs_ << " " << timestamp;
    fs_ << '\n';
    
    if ( P_.flush_records_ )
	fs_.flush();
  }

}

void nest::AsciiLogger2::record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  std::ofstream& fs_ = *(multimeters[id].fs_);
  #pragma omp critical
  {
    //print_id_(*multi_fs_, multi->multimeter_id);
    fs_ << id << ' ';
    fs_ << neuron_id << ' ';
    fs_ << timestamp;
    for (int i=0;i<data.size(); i++)
      fs_ << ' ' << data[i];
    fs_ << '\n';
    
    if ( P_.flush_records_ )
	fs_.flush();
    
  }

}

/*void nest::AsciiLogger2::append_value_to_multi_record(int id, int neuron_id, double v)
{
  #pragma omp critical
  {
    *multi_fs_ << ' ' << v;
    if (endrecord)
      *multi_fs_ << '\n';
  }
}*/

void nest::AsciiLogger2::signup_spike(int id, int neuron_id)
{
  #pragma omp critical
  {
    std::map<int,sd>::iterator it = spikedetectors.find(id);
    if(it != spikedetectors.end()) {
      it->second.neuron_id_.push_back(neuron_id);
    }
    else {
      spikedetectors.insert(std::pair<int,sd>(id, sd(id, new std::ofstream ())));
      spikedetectors[id].neuron_id_.push_back(neuron_id);
    }
    #ifdef _DEBUG_MODE
    std::cout << "signup_spike id=" << id << " neuron_id=" << neuron_id << std::endl;
    #endif
    //*spike_fs_ << spike->spikedetector_id << " " << neuron_id << "\n";
    //*spike_fs_ << "id=" << id << " neuron_id=" << " " << neuron_id << "\n";
  }
}

void nest::AsciiLogger2::signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames)
{
  #pragma omp critical
  {
    std::map<int,mm>::iterator it = multimeters.find(id);
    if(it != multimeters.end()) {
      it->second.neuron_id_.push_back(neuron_id);
    }
    else {
      multimeters.insert(std::pair<int,mm>(id, mm(id, new std::ofstream (), sampling_interval, valueNames)));
      multimeters[id].neuron_id_.push_back(neuron_id);
    }
    #ifdef _DEBUG_MODE
    std::cout << "signup_multi id=" << id << " neuron_id=" << neuron_id << " sampling_interval=" << sampling_interval << std::endl;
    #endif
  }
}

nest::AsciiLogger2::Parameters_::Parameters_(const std::string& path, const std::string& file_ext,
						bool withtime, bool withgid, bool withweight)
  : time_in_steps_(false),
    precise_times_(false),
    withgid_(withgid),
    withtime_(withtime),
    withweight_(withweight),
    precision_(3),
    scientific_(false),
    binary_(false),
    fbuffer_size_(BUFSIZ), // default buffer size as defined in <cstdio>
    label_(),
    file_ext_(file_ext),
    filename_(),
    close_after_simulate_(false),
    flush_after_simulate_(true),
    flush_records_(false),
    close_on_reset_(true),
    path_(path)
{}
nest::AsciiLogger2::AsciiLogger2(): n(0), P_(".", std::string("log"), true, true, false)
{  
  #ifdef _DEBUG_MODE
  std::cout << "nest::AsciiLogger2::AsciiLogger2()" << std::endl;
  #endif
}

nest::AsciiLogger2::AsciiLogger2(std::string path): n(0), P_(path, std::string("log"), true, true, false)
{  
}

void nest::AsciiLogger2::open_file(std::ofstream& fs_,int id, std::string prefix)
{  
  #ifdef _DEBUG_MODE
  std::cout << "nest::AsciiLogger2::open_file" << std::endl;
  #endif
  // we only close files here, opening is left to calibrate()
   if ( P_.close_on_reset_ && fs_.is_open() )
   {
     fs_.close();
     P_.filename_.clear();  // filename_ only visible while file open
   }
  
     // do we need to (re-)open the file
     bool newfile = false;

     if ( !fs_.is_open() )
     {
       newfile = true;   // no file from before
       P_.filename_ = build_filename_(id, prefix);
     }
     else
     {
       std::string newname = build_filename_(id, prefix);
       if ( newname != P_.filename_ )
       {
	  #ifdef NEST
         std::string msg = String::compose("Closing file '%1', opening file '%2'", P_.filename_, newname);
         Node::network()->message(SLIInterpreter::M_INFO, "RecordingDevice::calibrate()", msg);
         #endif

         fs_.close(); // close old file
         P_.filename_ = newname;
         newfile = true;
       }
     }

     if ( newfile )
     {
       #ifndef NESTIOPROXY
       assert(!fs_.is_open());
       #endif

       if ( P_.overwrite_files_ )
       {
         if ( P_.binary_ )
           fs_.open(P_.filename_.c_str(), std::ios::out | std::ios::binary);
         else
           fs_.open(P_.filename_.c_str());
       }
       else
       {
         // try opening for reading
         std::ifstream test(P_.filename_.c_str());
         if ( test.good() )
         {
	    #ifndef NESTIOPROXY
           std::string msg = String::compose("The device file '%1' exists already and will not be overwritten. "
                                             "Please change data_path, data_prefix or label, or set /overwrite_files "
                                             "to true in the root node.",P_.filename_);
           Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()", msg);
           throw IOError();
	    #endif
         }
         else
           test.close();

         // file does not exist, so we can open
         if ( P_.binary_ )
           fs_.open(P_.filename_.c_str(), std::ios::out | std::ios::binary);
         else
           fs_.open(P_.filename_.c_str());
       }

       if (P_.fbuffer_size_ != P_.fbuffer_size_old_)
       {
         if (P_.fbuffer_size_ == 0)
           fs_.rdbuf()->pubsetbuf(0, 0);
         else
         {
           std::vector<char>* buffer = new std::vector<char>(P_.fbuffer_size_);
           fs_.rdbuf()->pubsetbuf(reinterpret_cast<char*>(&buffer[0]), P_.fbuffer_size_);
         }
         
         P_.fbuffer_size_old_ = P_.fbuffer_size_;
       }
     }

     if ( !fs_.good() )
     {
       #ifndef NESTIOPROXY
       std::string msg = String::compose("I/O error while opening file '%1'. "
                                         "This may be caused by too many open files in networks "
					 "with many recording devices and threads.", P_.filename_);
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()", msg);
       #endif
       if ( fs_.is_open() )
         fs_.close();
       P_.filename_.clear();
       #ifndef NESTIOPROXY
       throw IOError();
       #endif
     }

     /* Set formatting
        Formatting is not applied to std::cout for screen output,
        since different devices may have different settings and
        this would lead to a mess.
      */
     if ( P_.scientific_ )
       fs_ << std::scientific;
     else
       fs_ << std::fixed;

     fs_ << std::setprecision(P_.precision_);

     
     #ifndef NESTIOPROXY
     if (P_.fbuffer_size_ != P_.fbuffer_size_old_)
     {
       std::string msg = String::compose("Cannot set file buffer size, as the file is already "
                                         "openeded with a buffer size of %1. Please close the "
                                         "file first.", P_.fbuffer_size_old_);
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()", msg);
       throw IOError();       
     }
     #endif
     n++;
}

void nest::AsciiLogger2::initialize(const double T) 
{
  #ifdef _DEBUG_MODE
  std::cout << "AsciiLogger2::initialize" << omp_get_thread_num() << std::endl;
  #endif
  int rank=0;
  //MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
  for (std::map<int,sd>::iterator it = spikedetectors.begin(); it!=spikedetectors.end(); ++it)
    open_file(*(it->second.fs_),it->first, "sd");
  
  
  for (std::map<int,mm>::iterator it = multimeters.begin(); it!=multimeters.end(); ++it)
    open_file(*(it->second.fs_),it->first, "mm");
  
  
  S_.files_open_ = true;
}

void nest::AsciiLogger2::close_file(std::ofstream& fs_)
{
  if ( fs_.is_open() )
   {
     if ( P_.close_after_simulate_ )
     {
       fs_.close();
       return;
     }

     if ( P_.flush_after_simulate_ )
       fs_.flush();

     if ( !fs_.good() )
     {
       #ifndef NESTIOPROXY
       std::string msg = String::compose("I/O error while opening file '%1'",P_.filename_);
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::finalize()", msg);

       throw IOError();
       #endif
     }
   }
}

void nest::AsciiLogger2::finalize()
{
  #ifdef _DEBUG_MODE
  std::cout << "AsciiLogger2::finalizeDatasets" << omp_get_thread_num() << std::endl;
  #endif
  for (std::map<int,sd>::iterator it = spikedetectors.begin(); it!=spikedetectors.end(); ++it)
    close_file(*(it->second.fs_));
  
  for (std::map<int,mm>::iterator it = multimeters.begin(); it!=multimeters.end(); ++it)
    close_file(*(it->second.fs_));
  
  S_.files_open_ = false;
}


nest::AsciiLogger2::~AsciiLogger2()
{
  if (S_.files_open_) {
    //print warning
    finalize();
  }
}

/*void AsciiLogger2::setBufferSize(int s)
{
    //
}*/

void nest::AsciiLogger2::syncronize(const double t)
{
    //
}

const std::string nest::AsciiLogger2::build_filename_(int id,std::string prefix) const
{
  // number of digits in number of virtual processes
#ifndef NESTIOPROXY
  
  Node& node = *(Node::network()->get_node(id));
  
  const int vpdigits = static_cast<int>(std::floor(std::log10(static_cast<float>(Communicator::get_num_virtual_processes()))) + 1);
  const int gidigits = static_cast<int>(std::floor(std::log10(static_cast<float>(Node::network()->size()))) + 1);
  std::ostringstream basename;
  const std::string& path = Node::network()->get_data_path();
  if (!P_.path_.empty() )
    basename << P_.path_ << '/';
  else if ( !path.empty() )
    basename << path << '/';
  basename << Node::network()->get_data_prefix();


  if ( !P_.label_.empty() )
    basename << P_.label_;
  else
    basename << node.get_name();

  basename << "-" << std::setfill('0') << std::setw(gidigits) << node.get_gid()
           << "-" << std::setfill('0') << std::setw(vpdigits) << node.get_vp(); //Zugriff nicht mehr möglich
  return basename.str() + '.' + P_.file_ext_;
#else
  const int vpdigits = omp_get_thread_num();
  const int gidigits = 0;
  
  std::ostringstream file_path;
  file_path << P_.path_ << '/';
  file_path << prefix << n << "." << P_.file_ext_;
  return file_path.str();
#endif
}

#ifndef NESTIOPROXY
void nest::AsciiLogger2::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;    // temporary copy in case of errors
  ptmp.set(S_,d);   // throws if BadProperty
  
  P_ = ptmp;
}
#endif


//copied from recording_device.cpp
//rd.mode_ lines are commented out
#ifndef NESTIOPROXY
void nest::AsciiLogger2::Parameters_::set(const State_& S, const DictionaryDatum& d)
{
  if (S.files_open_)
  {
    Node::network()->message(SLIInterpreter::M_WARNING, "AsciiLogger2::set_status", "Unpredictable behavior: Chaning ASCII Logger settings when files are open.");
  }
  
  updateValue<std::string>(d, names::label, label_);
  updateValue<bool>(d, names::withgid, withgid_);
  updateValue<bool>(d, names::withtime, withtime_);
  updateValue<bool>(d, names::withweight, withweight_);
  updateValue<bool>(d, names::time_in_steps, time_in_steps_);
  //if ( rd.mode_ == RecordingDevice::SPIKE_DETECTOR )
  //  updateValue<bool>(d, names::precise_times, precise_times_);
  updateValue<std::string>(d, names::file_extension, file_ext_);
  updateValue<long>(d, names::precision, precision_);
  updateValue<bool>(d, names::scientific, scientific_);

  updateValue<bool>(d, names::binary, binary_);

  long fbuffer_size;
  if (updateValue<long>(d, names::fbuffer_size, fbuffer_size))
  {  
    if (fbuffer_size < 0)
      throw BadProperty("/fbuffer_size must be <= 0");
    else
    {
      fbuffer_size_old_ = fbuffer_size_;
      fbuffer_size_ = fbuffer_size;
    }
  }
  
  updateValue<bool>(d, names::close_after_simulate, close_after_simulate_);
  updateValue<bool>(d, names::flush_after_simulate, flush_after_simulate_);
  updateValue<bool>(d, names::flush_records, flush_records_);
  updateValue<bool>(d, names::close_on_reset, close_on_reset_);

  // In Pynest we cannot use /record_to, because we have no way to pass
  // values as LiteralDatum. Thus, we must keep the boolean flags.
  // We must have || rec_change at the end, otherwise short-circuiting may
  // mean that some flags are not read.
}
#endif