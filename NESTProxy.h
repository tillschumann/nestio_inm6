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
#include <unistd.h>
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
		  if (s>0)
		    usleep(s);
		  sleeptimer[omp_get_thread_num()].pause();
		}
		
		void deliver(nestio::IDistribution* dist)
		{
		  delivertimer[omp_get_thread_num()].start();
		  int s = (int)dist->getValue();
		  //std::cout << "sleep " << s << std::endl;
		  if (s>0)
		    usleep(s);
		  delivertimer[omp_get_thread_num()].stop();
		}

		void sync(double& t)
		{
			logger.synchronize(t);
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
			  
			    nestio::IDistribution* numberOfSpikeDetectorsPerThread = conf.numberOfSpikeDetectorsPerThread->copy_init();
			    nestio::IDistribution* samplingIntervalsOfMeter = conf.samplingIntervalsOfMeter->copy_init();
			    nestio::IDistribution* numberOfValuesWrittenByMeter = conf.numberOfValuesWrittenByMeter->copy_init();
			    nestio::IDistribution* numberOfMultimetersPerThread = conf.numberOfMultimetersPerThread->copy_init();
			    
			    int thread_num = omp_get_thread_num();
			    bool spikedetector_config_for_thread = (conf.spikedetector_configs.find(nestio::getThreadHash())!=conf.spikedetector_configs.end());
			    int nosdpt;
			    if (spikedetector_config_for_thread)
			      nosdpt = conf.spikedetector_configs[nestio::getThreadHash()].size();
			    else
			      nosdpt = numberOfSpikeDetectorsPerThread->getIntValue();
			    
			    
			    
			    bool multimeter_config_for_thread = conf.multimeter_configs.find(nestio::getThreadHash())!=conf.multimeter_configs.end();
			    int nompt;
			    if (multimeter_config_for_thread) {
			      nompt = conf.multimeter_configs[nestio::getThreadHash()].size();
			    }
			    else {
			      nompt = numberOfMultimetersPerThread->getIntValue();
			    }
			    
			    #pragma omp critical 
			    {
			      spikeDetectors[thread_num].reserve(nosdpt);
			      multimeters[thread_num].reserve(nompt);
			    } 
			    
			    int neuron_id_offset = (thread_num*nosdpt)+rank*num_threads*nosdpt;
			    for (int i=0;i<nosdpt;i++) {
			      
				SpikeDetector* spikeDetector = new SpikeDetector(neuron_id_offset+i, &logger);
				if (spikedetector_config_for_thread) {
				  spikeDetector->connect2Neuron(neuron_id_offset+i,conf.spikedetector_configs[nestio::getThreadHash()].at(i).spikesPerDector->copy_init());
				}
				else {
				  spikeDetector->connect2Neuron(neuron_id_offset+i,conf.spikesPerDector->copy_init());
				}
				spikeDetector->singup();
				spikeDetectors[thread_num].push_back(spikeDetector);
			    }
			    
			    
			    int multimeter_id_offset = (thread_num*nompt)+rank*num_threads*nompt;
			    for (int i=0;i<nompt;i++) {
				double interval;
				Multimeter* multimeter;
				if (multimeter_config_for_thread) {
				  interval = conf.multimeter_configs[nestio::getThreadHash()].at(i).samplingIntervals->getValue();
				  multimeter = new Multimeter(multimeter_id_offset+i, interval, simSettings, conf.multimeter_configs[nestio::getThreadHash()].at(i).numberOfValuesWritten->getIntValue(), &logger);
				}
				else {
				  interval = samplingIntervalsOfMeter->getValue();
				  multimeter = new Multimeter(multimeter_id_offset+i, interval, simSettings, numberOfValuesWrittenByMeter->getIntValue(), &logger);

				}
				multimeter->connect2Neuron(multimeter_id_offset+i);
				multimeter->singup();
				multimeters[thread_num].push_back(multimeter);
			    }
			    
			    delete numberOfValuesWrittenByMeter;
			    delete samplingIntervalsOfMeter;
			    delete numberOfMultimetersPerThread;
			    delete numberOfSpikeDetectorsPerThread;
			    
			  
			  //#pragma omp barrier
			  
			}
			logger.initialize(simSettings.T);
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
			  nestio::IDistribution* deadTimeSpikeDetector = conf.deadTimeSpikeDetector->copy_init();
			  nestio::IDistribution* deadTimeMultimeters = conf.deadTimeMultimeters->copy_init();
			  nestio::IDistribution* deadTimeDeliver = conf.deadTimeDeliver->copy_init();
			
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
				    deliver(deadTimeDeliver);
				  }

				  {
				    // WRITE
				    //SCOREP_USER_REGION( "write", SCOREP_USER_REGION_TYPE_FUNCTION )
				    //std::cout << "update" << std::endl;
				    writetimer[thread_num].start();
				    for (int i=0;i<spikeDetectors[thread_num].size();i++) {
// 					//SCOREP_USER_REGION_BEGIN(epik_sleep,"sleep",SCOREP_USER_REGION_TYPE_COMMON);
					sleepAndMeasure(deadTimeSpikeDetector);
					//SCOREP_USER_REGION_END(epik_sleep);
					spikeDetectors[thread_num].at(i)->update(t,timestamp);
				    }
				    for (int i=0;i<multimeters[thread_num].size();i++) {
					//SCOREP_USER_REGION_BEGIN(epik_sleep,"sleep",SCOREP_USER_REGION_TYPE_COMMON);
					sleepAndMeasure(deadTimeMultimeters);
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
			
			delete deadTimeSpikeDetector;
			delete deadTimeMultimeters;
			delete deadTimeDeliver;
		    }
		    logger.finalize();
		}
};


#endif
