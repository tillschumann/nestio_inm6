#include "hdf5.h"
#include <cstring>
#include <string>
#include <vector>
#include "hdf5mpipp.h"
#include <random>

#ifndef NESTPROXY_CLASS
#define NESTPROXY_CLASS

class NESTProxy
{
	private:
		int iterations;
		void calc();
		void sync();
		void write();
		HDF5mpipp logger;
		
		int data;
	public:
		NESTProxy(int);
		~NESTProxy();
		void run();
};

class NESTProxyNode
{
	private:
		
	public:
		NESTProxyNode();
		~NESTProxyNode();
};

#endif