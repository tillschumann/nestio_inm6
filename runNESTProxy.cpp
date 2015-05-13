#include <mpi.h>

#include "NESTProxy.h"

#include "hdf5mpipp.h"
#include "ohdf5mpipp.h"
#include "sionlib_logger.h"
#include "AsciiLogger2.h"
//#include "abstract_logger.h"
#include "timer/scopetimer.h"
#include "timer/stopwatch.h"
#include "timer/seriestimer.h"

#include "nestio_func.h"

#include <iostream>
#include <sstream>
#include <omp.h>

#include <sys/types.h>
#include <sys/stat.h>

void run(nestio::Configuration &conf, nestio::SimSettings &simSettings,int argc, char *argv[])
{	
    int nop;
    MPI_Comm_size(MPI_COMM_WORLD,&nop);
    conf.numberOfProcesses = new nestio::FixIntValue(nop);
    
    
    std::stringstream output_dir_file;
    std::stringstream output_dir_log;
    if (argc > 2) {
      output_dir_file << argv[1];
      output_dir_log << argv[2];
    }
    else if (argc >1) {
      output_dir_file << argv[1];
      output_dir_log << argv[1];
    }
    else {
      //output_dir_file << "nestproxyoutput";
      //output_dir_log << "nestproxyoutput";
      
      std::cout << "Please set at least one argument for output folder location" << std::endl;
      std::cout << "if one argument is given: data and log files are written into this folder" << std::endl;
      std::cout << "if two arguemnts are give: first argument: data output folder, second argument: log output folder" << std::endl;
    }
    
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank==0) {
      mkdir(output_dir_file.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(output_dir_log.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    
    ILogger *logger;
    std::stringstream sfn, mfn;
    std::stringstream fns;
    
    switch (conf.logger) {
      case nestio::SIONLIB:
	//sfn << output_dir_file.str() << "/data.spike_sion";
	//mfn << output_dir_file.str() << "/data.multi_sion"; 
	logger = new Sionlib_logger(output_dir_file.str(), "sion", conf.bufferSize, conf.bufferSize, nestio::Standard);
	break;
      case nestio::SIONLIB_BUFFERED:
	//sfn << output_dir_file.str() << "/data.spike_sion";
	//mfn << output_dir_file.str() << "/data.multi_sion";
	logger = new Sionlib_logger(output_dir_file.str(), "sion", conf.bufferSize, conf.bufferSize, nestio::Buffered);
	break;
      case nestio::SIONLIB_COLLECTIVE:
	//sfn << output_dir_file.str() << "/data.spike_sion";
	//mfn << output_dir_file.str() << "/data.multi_sion";
	logger = new Sionlib_logger(output_dir_file.str(), "sion", conf.bufferSize, conf.bufferSize, nestio::Collective);
	break;
      case nestio::HDF5:
	fns << output_dir_file.str() << "/data.hdf5";
	logger = new HDF5mpipp(fns.str(), conf.bufferSize, simSettings);
	break;
      case nestio::oHDF5:
	fns << output_dir_file.str() << "/data_o.hdf5";
	logger = new OHDF5mpipp(fns.str(), conf.bufferSize, nestio::Standard);
	break;
      case nestio::oHDF5_BUFFERED:
	fns << output_dir_file.str() << "/data_o.hdf5";
	logger = new OHDF5mpipp(fns.str(), conf.bufferSize, nestio::Buffered);
	break;
      case nestio::oHDF5_COLLECTIVE:
	fns << output_dir_file.str() << "/data_o.hdf5";
	logger = new OHDF5mpipp(fns.str(), conf.bufferSize, nestio::Collective);
	break;
      case nestio::ASCII:
	sfn << output_dir_file.str();// << "/data.spike_ascii";
	mfn << output_dir_file.str();// << "/data.multi_ascii";
	logger = new nest::AsciiLogger2(sfn.str());
	break;
    }
    nest::SeriesTimer writetimer[conf.numberOfThreads->getIntValue()],
		      synctimer[conf.numberOfThreads->getIntValue()],
		      sleeptimer[conf.numberOfThreads->getIntValue()],
		      delivertimer[conf.numberOfThreads->getIntValue()];
    
    NESTProxy* proxy = new NESTProxy( simSettings,
		      conf,
		      *logger,
		      writetimer, synctimer, sleeptimer,delivertimer);
    proxy->run();
    
    std::stringstream benchmark_write_filename;
    std::stringstream benchmark_sync_filename;
    std::stringstream benchmark_sleep_filename;
    std::stringstream benchmark_deliver_filename;
    std::stringstream config_filename;
    
    benchmark_write_filename << output_dir_log.str() << "/benchfile_write.csv";
    benchmark_sync_filename << output_dir_log.str() << "/benchfile_sync.csv";
    benchmark_sleep_filename << output_dir_log.str() << "/benchfile_sleep.csv";
    benchmark_deliver_filename << output_dir_log.str() << "/benchfile_deliver.csv";
    config_filename << output_dir_log.str() << "/config.csv";
    
    
    if (rank==0) {
      std::ofstream configfile;
      configfile.open(config_filename.str(),std::ofstream::out|std::ofstream::app);
      configfile << conf.logger << ";";
      configfile << conf.bufferSize << ";";
      configfile << *(conf.numberOfProcesses) << ";";
      configfile << *(conf.numberOfThreads) << ";";
      configfile << *(conf.numberOfSpikeDetectorsPerThread) << ";";
      configfile << *(conf.spikesPerDector) << ";";
      configfile << *(conf.numberOfMultimetersPerThread) << ";";
      configfile << *(conf.samplingIntervalsOfMeter) << ";";
      configfile << *(conf.deadTimeSpikeDetector) << ";";
      configfile << *(conf.deadTimeMultimeters) << ";";
      configfile << *(conf.deadTimeDeliver) << ";";
      configfile << *(conf.numberOfValuesWrittenByMeter) << ";";
      configfile.close();
    }
    
    for(int r=0; r<conf.numberOfProcesses->getIntValue(); r++ ) {
      MPI_Barrier( MPI_COMM_WORLD );
      if (r==rank)
	for (int i=0;i<conf.numberOfThreads->getIntValue();i++) {
	  std::ofstream benchfile_write;
	  benchfile_write.open(benchmark_write_filename.str(),std::ofstream::out|std::ofstream::app);
	  std::stringstream ss2;
	  ss2 << rank << ";" << i;
	  writetimer[i].print_all_csv(ss2.str().c_str(),nest::Stopwatch::MILLISEC, benchfile_write);
	  benchfile_write.close();
	  
	  std::ofstream benchfile_sync;
	  benchfile_sync.open(benchmark_sync_filename.str(),std::ofstream::out|std::ofstream::app);
	  ss2.str("");
	  ss2 << rank << ";" << i;
	  synctimer[i].print_all_csv(ss2.str().c_str(), nest::Stopwatch::MILLISEC, benchfile_sync);
	  benchfile_sync.close();
	  
	  std::ofstream benchfile_sleep;
	  benchfile_sleep.open(benchmark_sleep_filename.str(),std::ofstream::out|std::ofstream::app);
	  ss2.str("");
	  ss2 << rank << ";" << i;
	  sleeptimer[i].print_all_csv(ss2.str().c_str(), nest::Stopwatch::MILLISEC, benchfile_sleep);
	  benchfile_sleep.close();
	  
	  std::ofstream benchfile_deliver;
	  benchfile_deliver.open(benchmark_deliver_filename.str(),std::ofstream::out|std::ofstream::app);
	  ss2.str("");
	  ss2 << rank << ";" << i;
	  delivertimer[i].print_all_csv(ss2.str().c_str(), nest::Stopwatch::MILLISEC, benchfile_deliver);
	  benchfile_deliver.close();
	}
    }
    delete proxy;
    delete logger;
}

int main(int argc, char *argv[])
{
        MPI_Init(&argc,&argv);
	
	#ifdef DEBUG_MODE
	std::cout << "runNESTProxy" << std::endl;
	#endif
	int numberOfThreads=omp_get_max_threads(); //must not be changed
	int threadNumber=omp_get_thread_num();
	
	nestio::SimSettings simSettings;
	simSettings.T=100.0;
	simSettings.Tresolution=0.1;
	
	nestio::Configuration conf;
	conf.logger = nestio::ASCII;
	conf.bufferSize = 2400;
	conf.numberOfThreads=new nestio::FixIntValue(numberOfThreads); //must not be changed
	conf.numberOfSpikeDetectorsPerThread=new nestio::FixIntValue(8);
	conf.spikesPerDector = new nestio::StandardDistribution(4.1,2.5);
	conf.numberOfMultimetersPerThread= new nestio::FixIntValue(8);
	conf.samplingIntervalsOfMeter = new nestio::FixDoubleValue(0.1);
	conf.deadTimeSpikeDetector = new nestio::FixIntValue(35); //
	conf.deadTimeMultimeters = new nestio::FixIntValue(35); //
	conf.deadTimeDeliver = new nestio::FixIntValue(100); //
	conf.numberOfValuesWrittenByMeter = new nestio::FixIntValue(4);
	
	
	//std::vector<nestio::Multimeter_Config> m_list_1_1(1);
	//m_list_1_1.at(0).numberOfValuesWritten=new nestio::FixIntValue(1);
	//m_list_1_1.at(0	).samplingIntervals=new nestio::FixDoubleValue(0.05);
	//std::pair<int,std::vector<nestio::Multimeter_Config>> pairvalue(nestio::getThreadHash(0,0),m_list_1_1);
	//conf.multimeter_configs.insert(pairvalue);

	run(conf,simSettings,argc,argv);
	#ifdef _DEBUG_MODE
	std::cout << "proxy complete" << std::endl;
	#endif
	delete conf.numberOfProcesses;
	delete conf.numberOfThreads;
	delete conf.numberOfSpikeDetectorsPerThread;
	delete conf.spikesPerDector;
	delete conf.numberOfMultimetersPerThread;
	delete conf.samplingIntervalsOfMeter;
	delete conf.deadTimeSpikeDetector;
	delete conf.deadTimeMultimeters;
	delete conf.deadTimeDeliver;
	delete conf.numberOfValuesWrittenByMeter;
	
	MPI_Finalize();
	
	
	#ifdef _DEBUG_MODE
	std::cout << "close" << std::endl;
	#endif
}
