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

#ifndef NESTPROXY_CLASS
#define NESTPROXY_CLASS

#define TWO_PI 6.2831853071795864769252866
 
/*double generateGaussianNoise(const double &variance)
{
	static bool hasSpare = false;
	static double rand1, rand2;
 
	if(hasSpare)
	{
		hasSpare = false;
		return sqrt(variance * rand1) * sin(rand2);
	}
 
	hasSpare = true;
 
	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
 
	return sqrt(variance * rand1) * cos(rand2);
}*/

template < typename L >
class NESTProxy
{
	private:
		L &logger;
		nest::SeriesTimer* writetimer;
		nest::SeriesTimer* synctimer;
		nest::SeriesTimer* sleeptimer;
		
		nestio::SimSettings simSettings;		
		nestio::Configuration conf;
		
		std::vector<SpikeDetector<L>*> spikeDetectors;
		std::vector<Multimeter<L>*> multimeters;
		
		
		void sleepAndMeasure(nestio::Distribution dist)
		{
		  sleeptimer[omp_get_thread_num()].start();
		  int s = (int)nestio::rand2(dist);
		  std::cout << "sleep " << s << std::endl;
		  sleep(s);
		  sleeptimer[omp_get_thread_num()].stop();
		}
		
		void deliver()
		{
		  sleepAndMeasure(conf.deadTimeDeliver);
		}

		/*void write(double& t)
		{
			int neuron_id=0;
			logger.single_write(t, data, neuron_id);
		}*/

		void sync()
		{
			logger.updateDatasetSizes();
		}
	public:
		NESTProxy(nestio::SimSettings &simSettings,
			  nestio::Configuration &conf,
			  L &logger,
			  nest::SeriesTimer* writetimer,
			  nest::SeriesTimer* synctimer,
			  nest::SeriesTimer* sleeptimer
 			):
		simSettings(simSettings),
		conf(conf),
		logger(logger),
		writetimer(writetimer),
		sleeptimer(sleeptimer),
		synctimer(synctimer)
		{
			std::cout << "Init NESTProxy" << std::endl;
			//set dataset name
			//only one dataset per node possible with hdf5
			int rank;
			MPI_Comm_rank (MPI_COMM_WORLD, &rank);
			
			//init rand
			sleep(rank);
			srand(time(NULL));
			std::cout << "print config" << std::endl;
			std::cout << "\tTstop="<< simSettings.T<< std::endl;
			std::cout << "\tTresolution="<< simSettings.Tresolution<< std::endl;
			std::cout << "configuration:" << std::endl;
			std::cout << conf << std::endl;
			

			//std::cout << "rank=" << rank << " conf.numberOfThreads=" << conf.numberOfThreads << " thread_num=" << thread_num << std::endl;
			int neuron_id_offset = (rank*conf.numberOfThreads)*conf.numberOfSpikeDetectorsPerThread;
			//std::cout << "neuron_id_offset=" << neuron_id_offset << std::endl;
			for (int i=0;i<conf.numberOfSpikeDetectorsPerThread;i++) {
			    SpikeDetector<L>* spikeDetector = new SpikeDetector<L>(neuron_id_offset+i, &logger);
			    spikeDetector->connect2Neuron(neuron_id_offset+i,conf.spikesPerDector);
			    spikeDetector->singup();
			    spikeDetectors.push_back(spikeDetector);
			}
			
			int multimeter_id_offset = (rank*conf.numberOfThreads)*conf.numberOfMultimetersPerThread;
			for (int i=0;i<conf.numberOfMultimetersPerThread;i++) {
			    double interval = nestio::rand2(conf.samlpingIntervalsOfMeter);
			    Multimeter<L>* multimeter = new Multimeter<L>(multimeter_id_offset+i, interval, simSettings, nestio::rand2(conf.numberOfValuesWrittenByMeter), &logger);
			    multimeter->connect2Neuron(multimeter_id_offset+i);
			    multimeter->singup();
			    multimeters.push_back(multimeter);
			}
			
			#pragma omp barrier
			logger.createDatasets();
		}

		~NESTProxy()
		{
		    std::cout << "destructor NESTProxy" << std::endl;
		    for (int i=0;i<conf.numberOfSpikeDetectorsPerThread;i++) {
			std::cout << "delete spikedetector " << i << std::endl;
			delete spikeDetectors.at(i);
		    }
		    for (int i=0;i<conf.numberOfMultimetersPerThread;i++) {
			std::cout << "delete multimeter " << i << std::endl;
			delete multimeters.at(i);
		    }
		    std::cout << "delete complete" << std::endl;
		}
		
		
		
		/*
		 * 
		 * Missing: a multimeter records more than one value per update step -> paramter mean/var mean + sqrt(var)
		 * 
		 */
		void run()
		{
			std::cout << "NEST PROXY RUN" << std::endl;
			std::cout << "Parameters:" << std::endl;
			std::cout << "RAND_MAX=" << RAND_MAX << std::endl;
			//std::cout << "thread_num=" << thread_num << std::endl;
			double t=simSettings.Tstart;
			int timestamp=0;
			
			//init
			//gather step
			
			
			while (t<=simSettings.T) {
				std::cout << "Iteration: t="<< t << " timestamp=" << timestamp << std::endl;
				// DELIVER
				//#pragma omp single
				
				#pragma omp parallel firstprivate(t,timestamp)
				{
				  int thread_num = omp_get_thread_num();
				  #pragma omp barrier
				  {
				    deliver();
				  }

				  {
				    // WRITE
				    writetimer[thread_num].start();
				    for (int i=0;i<conf.numberOfSpikeDetectorsPerThread;i++) {
					sleepAndMeasure(conf.deadTimeSpikeDetector);
					spikeDetectors.at(i)->update(t,timestamp);
				    }
				    for (int i=0;i<conf.numberOfMultimetersPerThread;i++) {
					sleepAndMeasure(conf.deadTimeMultimeters);
					multimeters.at(i)->update(t,timestamp);
				    }
				    writetimer[thread_num].stop();
				  }
				  //GATHER
				  #pragma omp barrier
				  //SYNC
				  synctimer[thread_num].start();
				  sync();
				  synctimer[thread_num].stop();
				
				}
				
				t+=simSettings.Tresolution;
				timestamp++;
			}
			
			std::cout << "Iterations done" << std::endl;
		}
};


#endif
