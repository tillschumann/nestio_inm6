#include <mpi.h>
#include <sion.h>
#include <vector>

#ifndef SIONLIBLOGGER_CLASS
#define SIONLIBLOGGER_CLASS

class Sionlib_logger
{
	private:
		int sid;
		FILE *fileptr;
		int own_id;
		sion_int64 buf_size;
		std::string prefix;
		std::vector<int> values;
		
	public:
		Sionlib_logger(std::string);
		~Sionlib_logger();
		void write(int data[]);
		void newDataSet(std::string, int);
		void setSize(int,int);
		void setBufferSize(int);
		void updateDatasetSizes();
		void single_write(int&);
};

#endif