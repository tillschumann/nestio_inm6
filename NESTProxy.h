#include "hdf5.h"
#include <cstring>
#include <string>
#include <vector>
#include "hdf5mpipp.h"
#include "sionlib_logger.h"
#include <random>
#include <mpi.h>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "timer/scopetimer.h"
#include "timer/stopwatch.h"
#include "timer/seriestimer.h"
#include "nestio_func.h"
#include <omp.h>

#include "Multimeter.h"
#include "SpikeDetector.h"

//#include <scorep/SCOREP_User.h>

#ifndef NESTPROXY_CLASS
#define NESTPROXY_CLASS


class NESTProxy
{
	private:
		ILogger &logger;
		nest::SeriesTimer* writetimer;
		nest::SeriesTimer* synctimer;
		nest::SeriesTimer* sleeptimer;
		nest::SeriesTimer* delivertimer;
		
		nestio::SimSettings simSettings;		
		nestio::Configuration conf;
		
		std::vector<SpikeDetector*>* spikeDetectors;
		std::vector<Multimeter*>* multimeters;
		
		int num_threads;
		
		
		void sleepAndMeasure(nestio::IDistribution* dist)
		{
		  sleeptimer[omp_get_thread_num()].start();
		  int s = (int)dist->getValue();
		  //std::cout << "sleep " << s << std::endl;
		  usleep(s);
		  sleeptimer[omp_get_thread_num()].pause();
		}
		
		void deliver()
		{
		  delivertimer[omp_get_thread_num()].start();
		  int s = (int)conf.deadTimeDeliver->getValue();
		  //std::cout << "sleep " << s << std::endl;
		  usleep(s);
		  delivertimer[omp_get_thread_num()].stop();
		}

		void sync(double& t)
		{
			logger.updateDatasetSizes(t);
		}
	public:
		NESTProxy(nestio::SimSettings &simSettings,
			  nestio::Configuration &conf,
			  ILogger &logger,
			  nest::SeriesTimer* writetimer,
			  nest::SeriesTimer* synctimer,
			  nest::SeriesTimer* sleeptimer,
			  nest::SeriesTimer* delivertimer
 			):
		simSettings(simSettings),
		conf(conf),
		logger(logger),
		delivertimer(delivertimer),
		writetimer(writetimer),
		sleeptimer(sleeptimer),
		synctimer(synctimer)
		{
			#ifdef _DEBUG_MODE
			std::cout << "Init NESTProxy" << std::endl;
			#endif
			//set dataset name
			//only one dataset per node possible with hdf5
			int rank;
			MPI_Comm_rank (MPI_COMM_WORLD, &rank);
			num_threads = omp_get_max_threads();
			
			//init rand
			sleep(rank);
			srand(time(NULL));
			#ifdef _DEBUG_MODE
			std::cout << "print config" << std::endl;
			std::cout << "\tTstop="<< simSettings.T<< std::endl;
			std::cout << "\tTresolution="<< simSettings.Tresolution<< std::endl;
			std::cout << "configuration:" << std::endl;
			std::cout << conf << std::endl;
			#endif
			
			
			
			spikeDetectors = new std::vector<SpikeDetector*>[num_threads];
			multimeters = new std::vector<Multimeter*>[num_threads];
			
			
			
			#pragma omp parallel 
			{
			    conf.numberOfSpikeDetectorsPerThread->init();
			    conf.numberOfMultimetersPerThread->init();
			    conf.samplingIntervalsOfMeter->init();
			    conf.numberOfValuesWrittenByMeter->init();
			    
			    int thread_num = omp_get_thread_num();
			    
			    int nosdpt = conf.numberOfSpikeDetectorsPerThread->getIntValue();
			    int nompt = conf.numberOfMultimetersPerThread->getIntValue();
			    
			    #pragma omp critical 
			    {
			      spikeDetectors[thread_num].reserve(nosdpt);
			      multimeters[thread_num].reserve(nompt);
			    }
			    
			    int neuron_id_offset = (thread_num*nosdpt)+rank*num_threads*nosdpt;
			    for (int i=0;i<nosdpt;i++) {
				SpikeDetector* spikeDetector = new SpikeDetector(neuron_id_offset+i, &logger);
				spikeDetector->connect2Neuron(neuron_id_offset+i,conf.spikesPerDector);
				spikeDetector->singup();
				spikeDetectors[thread_num].push_back(spikeDetector);
			    }
			    
			    
			    int multimeter_id_offset = (thread_num*nompt)+rank*num_threads*nompt;
			    for (int i=0;i<nompt;i++) {
				double interval = conf.samplingIntervalsOfMeter->getValue();
				Multimeter* multimeter = new Multimeter(multimeter_id_offset+i, interval, simSettings, conf.numberOfValuesWrittenByMeter->getIntValue(), &logger);
				multimeter->connect2Neuron(multimeter_id_offset+i);
				multimeter->singup();
				multimeters[thread_num].push_back(multimeter);
			    }
			  #pragma omp barrier
			  logger.createDatasets();
			}
		}

		~NESTProxy()
		{
		    #ifdef _DEBUG_MODE
		    std::cout << "destructor NESTProxy" << std::endl;
		    #endif
		    for (int thread_num=0;thread_num<num_threads;thread_num++) {
		      for (int i=0;i<spikeDetectors[thread_num].size();i++) {
			  #ifdef _DEBUG_MODE
			  std::cout << "delete spikedetector " << i << std::endl;
			  #endif
			  delete spikeDetectors[thread_num].at(i);
		      }
		      for (int i=0;i<multimeters[thread_num].size();i++) {
			  #ifdef _DEBUG_MODE
			  std::cout << "delete multimeter " << i << std::endl;
			  #endif
			  delete multimeters[thread_num].at(i);
		      }
		    }
		    #ifdef _DEBUG_MODE
		    std::cout << "delete complete" << std::endl;
		    #endif
		}
		
		
		
		/*
		 * 
		 * Missing: a multimeter records more than one value per update step -> paramter mean/var mean + sqrt(var)
		 * 
		 */
		void run()
		{
			#ifdef _DEBUG_MODE
			std::cout << "NEST PROXY RUN" << std::endl;
			std::cout << "Parameters:" << std::endl;
			std::cout << "RAND_MAX=" << RAND_MAX << std::endl;
			#endif
			//std::cout << "thread_num=" << thread_num << std::endl;
			double t=0;
			int timestamp=0;
			
			//init
			//gather step
			
			//SCOREP_USER_REGION_DEFINE(epik_sync);
			//SCOREP_USER_REGION_DEFINE(epik_sleep);
			
			#pragma omp parallel firstprivate(t,timestamp)
			{
			  conf.deadTimeSpikeDetector->init();
			  conf.deadTimeMultimeters->init();
			  conf.deadTimeDeliver->init();
			
			  while (t<=simSettings.T) {
				#ifdef _DEBUG_MODE
				std::cout << "Iteration: t="<< t << " timestamp=" << timestamp << std::endl;
				#endif
				// DELIVER
				//#pragma omp single
				  int thread_num = omp_get_thread_num();
				  #pragma omp barrier
				  {
				    //SCOREP_USER_REGION( "deliver", SCOREP_USER_REGION_TYPE_FUNCTION )
				    //std::cout << "deliver" << std::endl;
				    deliver();
				  }

				  {
				    // WRITE
				    //SCOREP_USER_REGION( "write", SCOREP_USER_REGION_TYPE_FUNCTION )
				    //std::cout << "update" << std::endl;
				    writetimer[thread_num].start();
				    for (int i=0;i<spikeDetectors[thread_num].size();i++) {
// 					//SCOREP_USER_REGION_BEGIN(epik_sleep,"sleep",SCOREP_USER_REGION_TYPE_COMMON);
					sleepAndMeasure(conf.deadTimeSpikeDetector);
					//SCOREP_USER_REGION_END(epik_sleep);
					spikeDetectors[thread_num].at(i)->update(t,timestamp);
				    }
				    for (int i=0;i<multimeters[thread_num].size();i++) {
					//SCOREP_USER_REGION_BEGIN(epik_sleep,"sleep",SCOREP_USER_REGION_TYPE_COMMON);
					sleepAndMeasure(conf.deadTimeMultimeters);
					//SCOREP_USER_REGION_END(epik_sleep);
					multimeters[thread_num].at(i)->update(t,timestamp);
				    }
				    //sleeptimer is paused and resumed in the sleepAndMeasure function
				    //to save the timings anyway stop is called
				    sleeptimer[thread_num].stop();
				    writetimer[thread_num].stop();
				  }
				  //GATHER
				  #pragma omp barrier
				  //SYNC
				  //SCOREP_USER_REGION_BEGIN(epik_sync,"sync",SCOREP_USER_REGION_TYPE_COMMON);
				  //std::cout << "sync" << std::endl;
				  synctimer[thread_num].start();
				  sync(t);
				  synctimer[thread_num].stop();
				  //SCOREP_USER_REGION_END(epik_sync);
				
				
				
				t+=simSettings.Tresolution;
				timestamp++;
				}
			#ifdef _DEBUG_MODE
			std::cout << "Iterations done" << std::endl;
			#endif
		  }
		}
};


#endif
