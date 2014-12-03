#include <mpi.h>

#include "NESTProxy.h"

#include "hdf5mpipp.h"
#include "sionlib_logger.h"
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

    ILogger *logger;
    switch (conf.logger) {
      case nestio::SIONLIB:
	logger = new Sionlib_logger("log.sion", conf.bufferSize, nestio::Standard, simSettings);
	break;
      case nestio::SIONLIB_BUFFERED:
	logger = new Sionlib_logger("log.sion", conf.bufferSize, nestio::Buffered, simSettings);
	break;
      case nestio::SIONLIB_COLLECTIVE:
	logger = new Sionlib_logger("log.sion", conf.bufferSize, nestio::Collective, simSettings);
	break;
      case nestio::HDF5:
	logger = new HDF5mpipp("log.hdf5", conf.bufferSize, simSettings);
	break;
    }
    nest::SeriesTimer writetimer[conf.numberOfThreads->getIntValue()],
		      synctimer[conf.numberOfThreads->getIntValue()],
		      sleeptimer[conf.numberOfThreads->getIntValue()],
		      delivertimer[conf.numberOfThreads->getIntValue()];
    
    NESTProxy proxy( simSettings,
		      conf,
		      *logger,
		      writetimer, synctimer, sleeptimer,delivertimer);
    proxy.run();
    
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    std::stringstream benchmark_write_filename;
    std::stringstream benchmark_sync_filename;
    std::stringstream benchmark_sleep_filename;
    std::stringstream benchmark_deliver_filename;
    std::stringstream config_filename;
    if (argc > 1) {
      benchmark_write_filename << "nestproxyoutput_" << argv[1];
      benchmark_sync_filename << "nestproxyoutput_" << argv[1];
      benchmark_sleep_filename << "nestproxyoutput_" << argv[1];
      benchmark_deliver_filename << "nestproxyoutput_" << argv[1];
      config_filename << "nestproxyoutput_" << argv[1];
    }
    else {
      benchmark_write_filename << "nestproxyoutput";
      benchmark_sync_filename << "nestproxyoutput";
      benchmark_sleep_filename << "nestproxyoutput";
      benchmark_deliver_filename << "nestproxyoutput";
      config_filename << "nestproxyoutput";
    }
    if (rank==0)
      mkdir(config_filename.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    benchmark_write_filename << "/benchfile_write.csv";
    benchmark_sync_filename << "/benchfile_sync.csv";
    benchmark_sleep_filename << "/benchfile_sleep.csv";
    benchmark_deliver_filename << "/benchfile_deliver.csv";
    config_filename << "/config.csv";
    
    
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
    delete logger;
}

int main(int argc, char *argv[])
{
        MPI_Init(&argc,&argv);
	
	
	std::cout << "runNESTProxy" << std::endl;
	int numberOfThreads=1;
	int threadNumber=0;
	#ifdef _OPENMP
	  numberOfThreads=omp_get_max_threads();
	  threadNumber = omp_get_thread_num();
	#endif
	
	
	nestio::SimSettings simSettings;
	simSettings.T=10.0;
	simSettings.Tresolution=0.1;
	
	nestio::Configuration conf;
	conf.logger = nestio::SIONLIB_BUFFERED;
	conf.bufferSize = 50;
	conf.numberOfThreads=new nestio::FixIntValue(numberOfThreads);
	conf.numberOfSpikeDetectorsPerThread=new nestio::FixIntValue(2);
	conf.spikesPerDector = new nestio::FixIntValue(1);
	conf.numberOfMultimetersPerThread= new nestio::FixIntValue(2);
	conf.samplingIntervalsOfMeter = new nestio::FixDoubleValue(0.05);
	conf.deadTimeSpikeDetector = new nestio::FixIntValue(35);
	conf.deadTimeMultimeters = new nestio::FixIntValue(35);
	conf.deadTimeDeliver = new nestio::FixIntValue(117426);
	conf.numberOfValuesWrittenByMeter = new nestio::FixIntValue(1);

	run(conf,simSettings,argc,argv);
	
	std::cout << "proxy complete" << std::endl;
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
	
	
	
	std::cout << "close" << std::endl;
}
