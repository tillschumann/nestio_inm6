#include <sion.h>
#include "iostream"
#include <cmath>
#include <map>
#include <vector>

#define HEADERSIZE 10

struct Multi  {
  int neuron_id;
  int numberOfValues;
  double interval;
  std::vector<std::string> valuesNames;
};

struct Boundaries {
 int min;
 int max;
};

bool inBoundaries(int value, Boundaries b)
{
 return (b.min <= value && b.max >= value);
}

Boundaries spikedetector_id_B {0,100};
Boundaries neuron_id_B {0,100};
Boundaries timestamp_B {0,5000};

Boundaries T_B{0,100};
Boundaries Tresolution_B{0,1};

Boundaries multimeter_id_B {0,100};
Boundaries numberOfValues_B{0,10};
Boundaries NodesCount_B{0,10};

Boundaries numberOfRecords_B{0,5000};
Boundaries values_B{0,100};


int readMultiFile(char* fn, bool printanyway)
{
  
  int errors=0;
  std::cout << "read sion file" << std::endl;
  
  //output varibales for sion_open function
  sion_int32 fsblksize;
  sion_int64 *chunksize = NULL;
  int *globalranks = NULL;
  int ntasks;
  int nfiles;
  FILE *fileptr;
  
  int sid = sion_open(fn, "rb", &ntasks, &nfiles, &chunksize, &fsblksize, &globalranks, &fileptr);
  
  //values are now know
  //std::cout << "ntasks=" << ntasks << std::endl;
  //std::cout << "fsblksize=" << fsblksize << std::endl;
  //std::cout << "chunksize[0]="<< chunksize[0] << std::endl;
  
  //iterating the nodes(tasks)
  for (int task=0; task<ntasks; task++) {
    int NodesCount;
    double Tstart;
    double T;
    int numberOfRecords=0;
    int startBody;
    double Tresolution;
    //seek to chunks of task
    sion_seek(sid, task, 0,0);
    
    //while (sion_feof(sid)<1)
    //  sion_fread(&numberOfRecords, sizeof(int), 1, sid);
    
    //sion_seek(sid, task, 0,0);
    
    //read values
    sion_fread(&NodesCount, sizeof(int), 1,sid);
    sion_fread(&T, sizeof(double), 1,sid);
    sion_fread(&Tresolution, sizeof(double), 1,sid);
    int numberOfRecords_wrong;
    sion_fread(&numberOfRecords_wrong, sizeof(int), 1,sid);
    sion_fread(&startBody, sizeof(int), 1,sid);
    
    
    if (printanyway) {
      std::cout << "NodesCount=" << NodesCount << std::endl;
      std::cout << "T=" << T << std::endl;
      std::cout << "Tresolution=" << Tresolution << std::endl;
      //std::cout << "numberOfRecords=" << numberOfRecords << std::endl;
      std::cout << "startBody=" << startBody << std::endl;
    }
    
    
    if (!inBoundaries(NodesCount, NodesCount_B)) {
      std::cerr << "NodesCount not in Boundaries" << std::endl;
      return errors+1;;
    }
    
    if (!inBoundaries(T, T_B)) {
      std::cerr << "T not in Boundaries" << std::endl;
      errors++;
    }
    
    if (!inBoundaries(Tresolution, Tresolution_B)) {
      std::cerr << "Tresolution not in Boundaries" << std::endl;
      errors++;
    }
    
    /*if (!inBoundaries(numberOfRecords, numberOfRecords_B)) {
      std::cerr << "numberOfRecords not in Boundaries" << std::endl;
      return errors+1;;
    }*/
    
    //read and store header information 
    std::map<int,Multi> idMap;
    for (int i=0; i<NodesCount; i++) {
      int multi_id;
      Multi multi;
     
      double interval;
      int numberOfValues;
      sion_fread(&multi_id, sizeof(int), 1,sid);
      sion_fread(&multi.neuron_id, sizeof(int), 1,sid);
      sion_fread(&multi.interval, sizeof(double), 1,sid);
      sion_fread(&multi.numberOfValues, sizeof(int), 1,sid);
      
      /*std::cout << "Node " << i << ":" << std::endl;
      std::cout << "\tmultimeter_id=" << multi_id << std::endl;
      std::cout << "\tneuron_id=" << multi.neuron_id << std::endl;
      std::cout << "\tinterval=" << multi.interval << std::endl;
      std::cout << "\tnumberOfValues=" << multi.numberOfValues << std::endl;*/
      
      if (!inBoundaries(multi.numberOfValues, numberOfValues_B)) {
	std::cerr << "multi.numberOfValues not in Boundaries" << std::endl;
	std::cerr << "multi.numberOfValues=" << multi.numberOfValues << std::endl;
	return errors+1;;
      }
      
      
      for (int v=0; v<multi.numberOfValues; v++) {
	char valueName[20];
	sion_fread(&valueName, 20, 1, sid);
	multi.valuesNames.push_back(valueName);
      }
      idMap.insert(std::pair<int,Multi>(multi_id,multi));

      /*std::cout << "\tvaluesNames=" << multi.valuesNames.at(0) ;
      for (int v=1; v<multi.numberOfValues; v++) {
	std::cout << "," << multi.valuesNames.at(v);
      }
      std::cout << std::endl;*/
    }
    
    double v[10];
    
    //read body of file
    
    //std::cout << "\t\tSTATUS\tM_ID\tN_ID\tTIMESTAMP\tVALUES" << std::endl;
    while (sion_feof(sid)<1) {
      int multi_id;
      sion_fread(&multi_id, sizeof(int), 1, sid);
      
      //if end of file multi_id contains numberOfRecords value
      if (sion_feof(sid)>0) {
	if (numberOfRecords != multi_id) {
	  std::cerr << "ERROR: numberOfRecords read != numberOfRecordss counted" << std::endl;
	  std::cerr << "numberOfRecords read ="<< multi_id <<std::endl;
	  std::cerr << "numberOfRecords counted ="<< numberOfRecords <<std::endl;
	  errors++;
	}
	break;
      }
      numberOfRecords++;
      int neuron_id;
      int numberOfValues;
      int timestamp;
      sion_fread(&neuron_id, sizeof(int), 1, sid);
      sion_fread(&timestamp, sizeof(int), 1, sid);
      sion_fread(&numberOfValues, sizeof(int), 1, sid);
      
      if (numberOfValues != idMap[multi_id].numberOfValues) {
	std::cerr << "ERROR: numberOfValues != idMap[multi_id].numberOfValues" << std::endl;
	std::cerr << "numberOfValues=" << numberOfValues << std::endl;
	std::cerr << "idMap[" << multi_id <<"].numberOfValues=" << idMap[multi_id].numberOfValues << std::endl;
      }
      sion_fread(&v, sizeof(double),idMap[neuron_id].numberOfValues,sid);
      
      bool valueNotInBoundaries=false;
      for (int i=1; i<idMap[multi_id].numberOfValues; i++) {
	  valueNotInBoundaries = (valueNotInBoundaries || inBoundaries(v[i], values_B));
      }
      
      if (!inBoundaries(multi_id, multimeter_id_B) ||
	 !inBoundaries(neuron_id, neuron_id_B) ||
	 !inBoundaries(numberOfValues, numberOfValues_B) ||
	 !inBoundaries(timestamp, timestamp_B) ||
	 valueNotInBoundaries || printanyway)
      {
	std::cerr << "OOB\t" << multi_id << "\t" << neuron_id << "\t" << timestamp << "\t";
	std::cerr << idMap[neuron_id].valuesNames.at(0) << "=" << v[0];
	for (int i=1; i<idMap[neuron_id].numberOfValues; i++) {
	  std::cerr << "\t" << idMap[neuron_id].valuesNames.at(i) << "=" << v[i];
	}
	std::cerr << std::endl;
	errors++;
      }
    }
  }

  
  //close file
  sion_close(sid);
  return errors;
}


int readSpikeFile(char* fn, bool printanyway)
{
  sion_int32 fsblksize;
  sion_int64 *chunksize = NULL;
  int *globalranks = NULL;
  int ntasks;
  int nfiles;
  int sid = sion_open(fn, "rb", &ntasks, &nfiles, &chunksize, &fsblksize, &globalranks, NULL);
  
  //std::cout << "ntasks: " << ntasks << std::endl;
  //std::cout << "chunksize[0]: " << chunksize[0] << std::endl;
  //std::cout << "fsblksize: " << fsblksize << std::endl;
  
  int errors=0;
  
  for (int task=0; task<ntasks; task++) {
    
    sion_seek(sid, task, 0,0);
    
    int nodesCount;
    double T;
    double Tresolution;
    int numberOfRecords=0;
    int startOfBody;
    
    //while (sion_feof(sid)<1)
    //  sion_fread(&numberOfRecords, sizeof(int), 1, sid);
    
    sion_seek(sid, task, 0,0);
    
    sion_fread(&nodesCount, sizeof(int), 1, sid);
    sion_fread(&T, sizeof(double), 1, sid);
    sion_fread(&Tresolution, sizeof(double), 1, sid); // should be Tresolution
    int numberOfRecords_wrong;
    sion_fread(&numberOfRecords_wrong, sizeof(int), 1, sid);
    sion_fread(&startOfBody, sizeof(int), 1, sid);
    
    /*std::cout << "task " << task << ":" << std::endl;
    std::cout << "\tchunksize: " << chunksize[task] << std::endl;
    std::cout << "\tglobalranks: " << globalranks[task] << std::endl;
    std::cout << std::endl;*/
    if (printanyway) {
      std::cout << "\tnodesCount: " << nodesCount << std::endl;
      std::cout << "\tT: " << T << std::endl;
      std::cout << "\tTresolution: " << Tresolution << std::endl;
      std::cout << "\tstartOfBody: " << startOfBody << std::endl;
      //std::cout << "\tnumberOfRecords: " << numberOfRecords << std::endl;
    }
    
    //int spikedetector_id;
    
    /*if (!inBoundaries(numberOfRecords, numberOfRecords_B)) {
      std::cerr << "numberOfRecords is not in Boundaries" << std::endl;
      std::cerr << "numberOfRecords=" << numberOfRecords << std::endl;
      return errors+1;
    }*/
    
    //std::cout << "\t\tSTATUS\tSD_ID\tN_ID\tTIMESTAMP" << std::endl;
    while (sion_feof(sid)<1) {
      int spikedetector_id;
      sion_fread(&spikedetector_id, sizeof(int), 1, sid);
      
      //if end of file multi_id contains numberOfRecords value
      if (sion_feof(sid)>0) {
	if (numberOfRecords != spikedetector_id) {
	  std::cerr << "ERROR: numberOfRecords read != numberOfRecordss counted" << std::endl;
	  std::cerr << "numberOfRecords read ="<< spikedetector_id <<std::endl;
	  std::cerr << "numberOfRecords counted ="<< numberOfRecords <<std::endl;
	  errors++;
	}
	break;
      }
      numberOfRecords++;
      
      int neuron_id;
      int timestamp;
      sion_fread(&neuron_id, sizeof(int), 1, sid);
      sion_fread(&timestamp, sizeof(int), 1, sid);
      
      if (!inBoundaries(spikedetector_id, spikedetector_id_B) ||
	 !inBoundaries(neuron_id, neuron_id_B) ||
	 !inBoundaries(timestamp, timestamp_B) || printanyway )	{
	 
	 std::cerr << "\t\tOOB\t";
	 std::cerr << spikedetector_id << "\t";
	 std::cerr << neuron_id << "\t";
	 std::cerr << timestamp << std::endl;
	
	 errors++;
	   
      }
    }
  }    
  
  sion_close(sid);
  return errors;
}

int main(int argc, char *argv[])
{
  int errors=0;
  bool printanyway=false;
  
  //look for options
  for (int i=0; i<argc; i++) {
    std::string argv_str(argv[i]);
    if (argv_str.find("printanyway")!=std::string::npos) {
      printanyway = true;
    }
  }
  
  //look for files
  for (int i=0; i<argc; i++) {
      std::string argv_str(argv[i]);
      
      if (argv_str.find(".spike_sion")!=std::string::npos) {
	std::cout << "spike file" << std::endl;
	errors += readSpikeFile(argv[i], printanyway);
      }
      
      if (argv_str.find(".multi_sion")!=std::string::npos) {
	std::cout << "multi file" << std::endl;
	errors += readMultiFile(argv[i], printanyway);
      }
  }
  
  if (errors>1 && !printanyway)
    std::cerr << errors << " errors found" << std::endl;
};