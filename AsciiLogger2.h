#ifndef ASCIILOGGER_CLASS
#define ASCIILOGGER_CLASS

#include <fstream>
//#include "multimeter.h"
//#include "spike_detector.h"
#include "abstract_logger.h"
#ifndef NESTIOPROXY
#include "recording_device.h"
#include "dictdatum.h"
#endif

/**
 * @file AsciiLogger2.h
 * Declarations for class AsciiLogger2.
 */

namespace nest
{
  class AsciiLogger2 : public ILogger
  {
	  private:
		  struct mm
		  {
		    const int id_;
		    mm(): id_(0) {};
		    mm(int id, std::ofstream *fs, double& sampling_interval, std::vector<Name>& valueNames): id_(id), fs_(fs), sampling_interval_(sampling_interval), valueNames_(valueNames) {};
		       
		    std::vector<int> neuron_id_;
		    double sampling_interval_;
		    std::vector<Name> valueNames_;
		    
		    std::ofstream *fs_;
		    
		    bool operator<( const mm & n ) const {
		      return this->id_ < n.id_;
		    }
		  };
		  
		  
		  
		  struct sd
		  {
		    const int id_;
		    sd(): id_(0) {};
		    sd(int id, std::ofstream *fs): id_(id), fs_(fs) {};
		    
		    std::vector<int> neuron_id_;
		    std::ofstream *fs_;
		    
		    
		    bool operator<( const sd & n ) const {
		      return this->id_ < n.id_;
		    }
		  };
		  
		  
		  std::map<int,sd> spikedetectors;
		  std::map<int,mm> multimeters;

		  //std::ofstream spike_fs_;
		  //std::ofstream multi_fs_;
		  
		  //counts open file times
		  int n;
		  
		  struct State_ {
		    bool files_open_;
		    
		    State_(): files_open_(false) {};
		  };
		  
		  struct Parameters_ {
		    bool time_in_steps_; //!< true if time is printed in steps, not ms.
		    bool precise_times_; //!< true if time is computed including offset
		    bool withgid_;       //!< true if element GID is to be printed, default
		    bool withtime_;      //!< true if time of event is to be printed, default
		    bool withweight_;    //!< true if weight of event is to be printed

		    long precision_;     //!< precision of doubles written to file
		    bool scientific_;    //!< use scientific format if true, else fixed

		    bool binary_;            //!< true if to write files in binary mode instead of ASCII   
		    long fbuffer_size_;      //!< the buffer size to use when writing to file
		    long fbuffer_size_old_;  //!< the buffer size to use when writing to file (old)

		    std::string label_;    //!< a user-defined label for symbolic device names.
		    std::string file_ext_; //!< the file name extension to use, without .
		    std::string filename_; //!< the filename, if recording to a file (read-only)
		    bool close_after_simulate_;    //!< if true, finalize() shall close the stream
		    bool flush_after_simulate_;    //!< if true, finalize() shall flush the stream
		    bool flush_records_;   //!< if true, flush stream after each output
		    bool close_on_reset_;  //!< if true, close stream in init_buffers() 
		    
		    
		    bool overwrite_files_;
		    std::string path_;

		    /**
		    * Set default parameter values.
		    * @param Default file name extension, excluding ".".
		    * @param Default value for withtime property
		    * @param Default value for withgid property
		    */
		    Parameters_(const std::string&, const std::string&, bool, bool, bool);

		    //void get(const AsciiLogger2&, DictionaryDatum&) const;  //!< Store current values in dictionary
		    #ifndef NESTIOPROXY
		    void set(const State_&, const DictionaryDatum&);  //!< Set values from dicitonary
		    #endif
		  };
		 
		  
		  Parameters_ P_;
		  State_ S_;
		  
		   /**
		    * Build filename from parts.
		    * @note This function returns the filename, it does not manipulate
		    *       any data member.
		    */
		    const std::string build_filename_(int id,std::string prefix) const;
		    
		    void open_file(std::ofstream& fs_, int id,  std::string prefix);
		    void close_file(std::ofstream& fs_);
		  
		  
	  public:
		  AsciiLogger2();
		  AsciiLogger2(std::string path);
		  ~AsciiLogger2();
		  //int newDataSet(const std::string, const int, const int);
		  //void setSize(int,int);
		  //void setBufferSize(int);
		  //void single_write(double& t, int& v, const int ptr);
		  //void single_write(double& t, double& v, const int ptr);
		  
		  
		  //void get_status(DictionaryDatum &) const; TODO-IO
		  #ifndef NESTIOPROXY
		  void set_status(const DictionaryDatum &);
		  #endif
		  
		  void record_spike(int id, int neuron_id, int timestamp);
		  void record_multi(int id, int neuron_id, int timestamp, const std::vector<double_t>& data);
		  //void append_value_to_multi_record(int id, int neuron_id, double v, bool endrecord);
		  void signup_spike(int id, int neuron_id);
		  
		  void signup_multi(int id, int neuron_id, double sampling_interval, std::vector<Name> valueNames);
		  
		  
		  void syncronize(const double t);
		  void initialize(const double T);
		  void finalize();
  };
}

#endif
