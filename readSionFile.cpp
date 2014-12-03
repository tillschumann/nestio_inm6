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


int readMultiFile()
{
  std::cout << "read sion file" << std::endl;
  
  //output varibales for sion_open function
  sion_int32 fsblksize;
  sion_int64 *chunksize = NULL;
  int *globalranks = NULL;
  int ntasks;
  int nfiles;
  FILE *fileptr;
  
  int sid = sion_open("multi_log.sion", "rb", &ntasks, &nfiles, &chunksize, &fsblksize, &globalranks, &fileptr);
  
  //values are now know
  std::cout << "ntasks=" << ntasks << std::endl;
  std::cout << "fsblksize=" << fsblksize << std::endl;
  std::cout << "chunksize[0]="<< chunksize[0] << std::endl;
  
  int NodesCount;
  double Tstart;
  double T;
  int numberOfRecords;
  int startBody;
  double Tresolution;
  
  //iterating the nodes(tasks)
  for (int task=0; task<ntasks; task++) {
    
    //seek to chunks of task
    sion_seek(sid, task, 0,0);
    
    //read values
    sion_fread(&NodesCount, sizeof(int), 1,sid);
    sion_fread(&T, sizeof(double), 1,sid);
    sion_fread(&Tresolution, sizeof(double), 1,sid);
    sion_fread(&numberOfRecords, sizeof(int), 1,sid);
    sion_fread(&startBody, sizeof(int), 1,sid);
    
    std::cout << "NodesCount=" << NodesCount << std::endl;
    std::cout << "T=" << T << std::endl;
    std::cout << "Tresolution=" << Tresolution << std::endl;
    std::cout << "numberOfRecords=" << numberOfRecords << std::endl;
    std::cout << "startBody=" << startBody << std::endl;
    
    
    //read and store header information 
    std::map<int,Multi> idMap;
    for (int i=0; i<NodesCount; i++) {
      int multi_id;
      Multi multi;
     
      double interval;
      int numberOfValues;
      sion_fread(&multi_id, sizeof(int), 1,sid);
      sion_fread(&multi.interval, sizeof(double), 1,sid);
      sion_fread(&multi.numberOfValues, sizeof(int), 1,sid);
      for (int v=0; v<multi.numberOfValues; v++) {
	char valueName[20];
	sion_fread(&valueName, 20, 1, sid);
	multi.valuesNames.push_back(valueName);
      }
      idMap.insert(std::pair<int,Multi>(multi_id,multi));
      
      std::cout << "Node " << i << ":" << std::endl;
      std::cout << "\tmultimeter_id=" << multi_id << std::endl;
      //std::cout << "\tneuron_id=" << multi.neuron_id << std::endl;
      std::cout << "\tinterval=" << multi.interval << std::endl;
      std::cout << "\tnumberOfValues=" << multi.numberOfValues << std::endl;
      std::cout << "\tvaluesNames=" << multi.valuesNames.at(0) ;
      for (int v=1; v<multi.numberOfValues; v++) {
	std::cout << "," << multi.valuesNames.at(v);
      }
      std::cout << std::endl;
    }
    
    double v[10];
    std::cout << "BODY" << std::endl;
    
    //read body of file
    while (sion_feof(sid)<1) {
      int multi_id;
      int neuron_id;
      int numberOfValues;
      int timestamp;
      sion_fread(&multi_id, sizeof(int), 1, sid);
      sion_fread(&neuron_id, sizeof(int), 1, sid);
      sion_fread(&timestamp, sizeof(int), 1, sid);
      sion_fread(&numberOfValues, sizeof(int), 1, sid);
      
      if (numberOfValues != idMap[multi_id].numberOfValues)
	std::cout << "ERROR: numberOfValues != idMap[multi_id].numberOfValues" << std::endl;
      
      std::cout << "multimeter_id=" << multi_id << "\tneuron_id=" << neuron_id << "\ttimestamp=" << timestamp << "\t";
      sion_fread(&v, sizeof(double),idMap[multi_id].numberOfValues,sid);
      std::cout << idMap[multi_id].valuesNames.at(0) << "=" << v[0];
      for (int i=1; i<idMap[multi_id].numberOfValues; i++) {
	std::cout << "\t" << idMap[multi_id].valuesNames.at(i) << "=" << v[i];
      }
      std::cout << std::endl;
    }
  }

  
  //close file
  sion_close(sid);
  return 0;
}


int readSpikeFile()
{
  sion_int32 fsblksize;
  sion_int64 *chunksize = NULL;
  int *globalranks = NULL;
  int ntasks;
  int nfiles;
  int sid = sion_open("spikes_log.sion", "rb", &ntasks, &nfiles, &chunksize, &fsblksize, &globalranks, NULL);
  
  std::cout << ntasks << std::endl;
  
  for (int task=0; task<ntasks; task++) {
    
    int nodesCount;
    double T;
    double Tresolution;
    int numberOfRecords;
    int startOfBody;
    
    sion_fread(&nodesCount, sizeof(int), 1, sid);
    sion_fread(&simSettings.T, sizeof(double), 1, sid);
    sion_fread(&simSettings.Tresolution, sizeof(double), 1, sid); // should be Tresolution
    sion_fread(&numberOfRecords, sizeof(int), 1, sid);
    sion_fread(&startOfBody, sizeof(int), 1, sid);
    
    
    
  }
  
  
  sion_close(sid);
  return 0;
}

int main(int argc, char *argv[])
{
  readSpikeFile();
};