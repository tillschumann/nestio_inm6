#include <sion.h>
#include "iostream"
#include <cmath>

#define HEADERSIZE 10

int main(int argc, char *argv[])
{
  std::cout << "read sion file" << std::endl;
  sion_int32 fsblksize;
  sion_int64 *chunksize = NULL;
  int *globalranks = NULL;
  int ntasks;
  int nfiles;
  FILE *fileptr;
  
  int sid = sion_open("multi_log.sion", "rb", &ntasks, &nfiles, &chunksize, &fsblksize, &globalranks, &fileptr);
  
  std::cout << "ntasks=" << ntasks << std::endl;
  std::cout << "fsblksize=" << fsblksize << std::endl;
  
  int NodesCount;
  double Tstart;
  double T;
  
  sion_fread(&NodesCount, sizeof(int), 1,sid);
  sion_fread(&Tstart, sizeof(double), 1,sid);
  sion_fread(&T, sizeof(double), 1,sid);
  
  std::cout << "NodesCount=" << NodesCount << std::endl;
  std::cout << "Tstart=" << Tstart << std::endl;
  std::cout << "T=" << T << std::endl;
  
  for (int t=0; t<ntasks; t++) {
    
    sion_seek(sid, t, 0, 0); // only for nfiles <= 1
    for (int i=0; i<NodesCount; i++) {
      int id;
      int owner_id;
      double interval;
      sion_fread(&id, sizeof(int), 1,sid);
      sion_fread(&owner_id, sizeof(int), 1,sid);
      sion_fread(&interval, sizeof(double), 1,sid);
      
      std::cout << "Node " << i << ":" << std::endl;
      std::cout << "\tid=" << id << std::endl;
      std::cout << "\towner_id=" << owner_id << std::endl;
      std::cout << "\tinterval=" << interval << std::endl;
    }
    int chunknum = (int)ceil(8*HEADERSIZE/ chunksize[t]);
    int posinchunk = 8*HEADERSIZE% chunksize[t];
    sion_seek(sid, t, chunknum, posinchunk); // only for nfiles <= 1
  }
  while ()
  
  
  
  sion_close(sid);
};