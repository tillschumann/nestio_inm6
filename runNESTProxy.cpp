#include <mpi.h>

#include "NESTProxy.h"

#include "hdf5mpipp.h"
#include "sionlib_logger.h"
#include "timer/scopetimer.h"
#include "timer/stopwatch.h"
#include "timer/seriestimer.h"

#include "nestio_func.h"

#include <iostream>
#include <sstream>
#include <omp.h>

int main(int argc, char *argv[])
{
	std::cout << "runNESTProxy" << std::endl;
	int numberOfThreads=1;
	/*#pragma omp parallel 
	{
	  numberOfThreads=omp_get_num_threads();
	  std::cout << "number of threads: " << numberOfThreads << std::endl;
	}
	#endif
	*/
	int threadNumber=0;
	#ifdef _OPENMP
	  numberOfThreads=omp_get_max_threads();
	  threadNumber = omp_get_thread_num();
	#endif
	
	
	nestio::SimSettings simSettings;
	
	simSettings.Tstart=0;
	simSettings.T=0.5;
	simSettings.Tresolution=0.1;
	
	nestio::Configuration conf;
	conf.numberOfThreads=numberOfThreads;
	conf.numberOfProcesses=1;
	conf.numberOfSpikeDetectorsPerThread=3;
	conf.spikesPerDector.mean=5;
	conf.spikesPerDector.var=3;
	conf.numberOfMultimetersPerThread=2;
	conf.samlpingIntervalsOfMeter.mean=0.2;
	conf.samlpingIntervalsOfMeter.var=0;
	conf.deadTimeSpikeDetector.mean=0;
	conf.deadTimeSpikeDetector.var=0;
	conf.deadTimeMultimeters.mean=0;
	conf.deadTimeMultimeters.var=0;
	conf.deadTimeDeliver.mean=0;
	conf.deadTimeDeliver.var=0;
	conf.numberOfValuesWrittenByMeter.mean=3;
	conf.numberOfValuesWrittenByMeter.var=2;
	
	MPI_Init(&argc,&argv);
	{	
	  
		MPI_Comm_size(MPI_COMM_WORLD,&conf.numberOfProcesses);
		//HDF5mpipp logger_hdf5("log.hdf5", 1000, simSettings);
		Sionlib_logger logger_sion("log.sion", 1000, simSettings);
		
		  
		nest::SeriesTimer writetimer[conf.numberOfThreads], synctimer[conf.numberOfThreads], sleeptimer[conf.numberOfThreads];

		
		/*NESTProxy<HDF5mpipp> proxy(simSettings,
					    conf,
					   logger_hdf5,
					   writetimer, synctimer, sleeptimer);
		proxy.run();*/
		
		
		NESTProxy<Sionlib_logger> proxy2(simSettings,
						 conf,
						 logger_sion,
						 writetimer, synctimer, sleeptimer);
		proxy2.run();
		
		int rank;
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
		for(int r=0; r<conf.numberOfProcesses; r++ ) {
		  MPI_Barrier( MPI_COMM_WORLD );
		  if (r==rank)
		    for (int i=0;i<conf.numberOfThreads;i++) {
		      std::ofstream benchfile;
		      //std::stringstream ss;
		      //ss << "benchfile_" << rank << "_" << i;
		      benchfile.open("benchfile.csv",std::ofstream::out|std::ofstream::app);
		      
		      std::stringstream ss2;
		      ss2 << rank << ";" << i  << ";write timings";
		      writetimer[i].print_csv(ss2.str().c_str(),nest::Stopwatch::MILLISEC, benchfile);
		      ss2.str("");
		      ss2 << rank << ";" << i << ";sync timings";
		      synctimer[i].print_csv(ss2.str().c_str(), nest::Stopwatch::MILLISEC, benchfile);
		      ss2.str("");
		      ss2 << rank << ";" << i << ";sleep timings" ;
		      sleeptimer[i].print_csv(ss2.str().c_str(), nest::Stopwatch::MILLISEC, benchfile);
		      benchfile.close();
		    }
		}
	}
	
	std::cout << "proxy complete" << std::endl;

	
	MPI_Finalize();
	
	std::cout << "close" << std::endl;
}
