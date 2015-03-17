#include "AsciiLogger2.h"
#include "iostream"
#include "sstream"
#include <fstream> 
#include <iomanip>
//#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

void nest::AsciiLogger2::record_spike(int id, int neuron_id, int timestamp)
{
  for (int i=0; i<spikedetectors.size();i++) {
    if (spikedetectors[i].id == id && spikedetectors[i].neuron_id == neuron_id) {
      std::ofstream& fs_ = *(spikedetectors[i].fs_);
      #pragma omp critical
      {
	fs_ << id;
	fs_ << " " << neuron_id;
	fs_ << " " << timestamp;
	fs_ << '\n';
	
	if ( P_.flush_records_ )
	    fs_.flush();
      }
      break;
    }
  }
}

void nest::AsciiLogger2::record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data)
{
  for (int i=0; i<multimeters.size();i++) {
    if (multimeters[i].id == id && multimeters[i].neuron_id == neuron_id) {
      std::ofstream& fs_ = *(spikedetectors[i].fs_);
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
      break;
    }
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

void nest::AsciiLogger2::signup_spike(int id, int neuron_id, int expectedSpikeCount)
{
  #pragma omp critical
  {
    spikedetectors.push_back((sd){id, neuron_id, new std::ofstream ()});
    std::cout << "signup_spike id=" << id << " neuron_id=" << neuron_id << " expectedSpikeCount=" << expectedSpikeCount << std::endl;
    //*spike_fs_ << spike->spikedetector_id << " " << neuron_id << "\n";
    //*spike_fs_ << "id=" << id << " neuron_id=" << " " << neuron_id << "\n";
  }
}

void nest::AsciiLogger2::signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames)
{
  #pragma omp critical
  {
    multimeters.push_back((mm){id, neuron_id, new std::ofstream ()});  //TODO
    std::cout << "signup_multi id=" << id << " neuron_id=" << neuron_id << " sampling_interval=" << sampling_interval << std::endl;
    //*multi_fs_ << "id=" << id << " neuron_id=" << " " << neuron_id << " " << "-1" << " ";
    //*multi_fs_ << multi->multimeter_id << " " << neuron_id << " " << multi->numberOfValues << " ";
    //for (int i=0; i<multi->numberOfValues; i++)
    //  *multi_fs_ << multi->valueNames[i] << " ";
    //*multi_fs_ << "\n";
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

nest::AsciiLogger2::AsciiLogger2(std::string path): n(0), P_(path, std::string("log"), true, true, false)
{  
}

void nest::AsciiLogger2::open_file(std::ofstream& fs_, std::string prefix)
{  
  std::cout << "nest::AsciiLogger2::open_file" << std::endl;
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
       P_.filename_ = build_filename_(prefix);
     }
     else
     {
       std::string newname = build_filename_(prefix);
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
       #ifdef NEST
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
	    #ifdef NEST
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
       #ifdef NEST
       std::string msg = String::compose("I/O error while opening file '%1'. "
                                         "This may be caused by too many open files in networks "
					 "with many recording devices and threads.", P_.filename_);
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()", msg);
       #endif
       if ( fs_.is_open() )
         fs_.close();
       P_.filename_.clear();
       #ifdef NEST
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

     
     #ifdef NEST
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
  std::cout << "AsciiLogger2::createDatasets" << omp_get_thread_num() << std::endl;
  std::cout << "T=" << T << std::endl;
  int rank=0;
  //MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  
  for (int i=0;i<spikedetectors.size();i++)
    open_file(*(spikedetectors[i].fs_), "sd");
  
  for (int i=0;i<multimeters.size();i++)
    open_file(*(multimeters[i].fs_), "mm");
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
       #ifdef NEST
       std::string msg = String::compose("I/O error while opening file '%1'",P_.filename_);
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::finalize()", msg);

       throw IOError();
       #endif
     }
   }
}

void nest::AsciiLogger2::finalize()
{
  std::cout << "AsciiLogger2::finalizeDatasets" << omp_get_thread_num() << std::endl;
  for (int i=0;i<spikedetectors.size();i++)
    close_file(*(spikedetectors[i].fs_));
  
  for (int i=0;i<multimeters.size();i++)
    close_file(*(multimeters[i].fs_));
}


nest::AsciiLogger2::~AsciiLogger2()
{
    
}

/*void AsciiLogger2::setBufferSize(int s)
{
    //
}*/

void nest::AsciiLogger2::syncronize(const double t)
{
    //
}

const std::string nest::AsciiLogger2::build_filename_(std::string prefix) const
{
  // number of digits in number of virtual processes
#ifdef NEST
  const int vpdigits = static_cast<int>(std::floor(std::log10(static_cast<float>(Communicator::get_num_virtual_processes()))) + 1);
  const int gidigits = static_cast<int>(std::floor(std::log10(static_cast<float>(Node::network()->size()))) + 1);
  std::ostringstream basename;
  const std::string& path = Node::network()->get_data_path();
  if ( !path.empty() )
    basename << path << '/';
  basename << Node::network()->get_data_prefix();


  if ( !P_.label_.empty() )
    basename << P_.label_;
  else
    basename << node_.get_name();

  basename << "-" << std::setfill('0') << std::setw(gidigits) << node_.get_gid()
           << "-" << std::setfill('0') << std::setw(vpdigits) << node_.get_vp();
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
