#include "NESTProxy.h"

#include <mpi.h>
#include <sstream>
#include <random>
#include <math.h>
#include <stdlib.h>

NESTProxy::NESTProxy(int iterations): iterations(iterations), logger("log.h5")
{
	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	std::stringstream ss;
	ss << "N" << rank;
    logger.newDataSet(ss.str(), 1000);
	logger.setBufferSize(1);
}

NESTProxy::~NESTProxy()
{
}

void NESTProxy::run()
{
	int next=0;
    int n=0;
	while (n<iterations) {
		sync();
		calc();
		if (n>=next) {
			next += logf(1.0f - (float) random() / (RAND_MAX + 1)) * (iterations/10.);
			std::cout << "next=" << next << std::endl;
			write();
		}
		n++;
	}
}

void NESTProxy::calc()
{
  int x;
  MPI_Comm_rank (MPI_COMM_WORLD, &x);
  data += x;
}

void NESTProxy::write()
{
	std::cout << "write data" << std::endl;
	logger.write(&data);
}

void NESTProxy::sync()
{
	std::cout << "sync" << std::endl;
	logger.updateDatasetSizes();
}