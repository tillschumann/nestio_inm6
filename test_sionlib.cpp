#include <mpi.h>
#include <sion.h>


int main(int argc, char *argv[])
{
	/* SION parameters */
	int sid, numFiles, globalrank, my_rank, num_procs;
	MPI_Comm lComm;
	sion_int64 chunksize, left, bwrote;
	sion_int32 fsblksize;
	char fname[256], *newfname=NULL;
	FILE *fileptr;
	/* initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	/* open parameters */
	chunksize = 10; globalrank = my_rank;
	strcpy(fname, "parfile.sion");
	numFiles = 1; fsblksize = -1;
	
	/* create a new file */
	sid = sion_paropen_mpi(fname, "bw", &numFiles,
	MPI_COMM_WORLD, &lComm,
	&chunksize, &fsblksize,
	&globalrank,
	&fileptr, &newfname);
	/* write buffer to file */
	left=chunksize;
	char* p = (char *) fname;
	while (left > 0) {
		sion_ensure_free_space(sid, left);
		bwrote = fwrite(p, 1, left, fileptr);
		left -= bwrote; p += bwrote;
	}
	/* close file */
	sion_parclose_mpi(sid);
	/* finalize MPI */
	MPI_Finalize();
}