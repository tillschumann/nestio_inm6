#include <mpi.h>

#include "NESTProxy.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	{
		NESTProxy proxy(100);
		proxy.run();
	}
	MPI_Finalize();
}