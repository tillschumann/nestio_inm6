#include "nestio_func.h"

int main( int argc, const char* argv[] )
{
	nestio::Distribution dist(187881766,85435260);
	std::cout << "dist:" << std::endl;
	double min, max;
	min=187881766;
	max=min;
	for (int i=0; i<100; i++) {
	  double v = dist.getValue();
	  std::cout << v << std::endl;
	  if (v<min) min=v;
	  if (v>max) max=v;
	}
	
	std::cout << "MAX=" << max << std::endl;
	std::cout << "MIN=" << min << std::endl;
}