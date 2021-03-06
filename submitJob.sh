#!/usr/bin/env bash

# Remarks: A line beginning with # is a comment.
#          A line beginning with #PBS is a PBS directive.
#          PBS directives must come first; any directives
#             after the first executable statement are ignored.
#

### The path and file name for standard output.
#PBS -o /users/schumann/test/nestio/job_out

### The path and file name for error output.
#PBS -e /users/schumann/test/nestio/job_out

### set name of job
#PBS -N hdf5experiment

### set the number of nodes and processes per node (ppn)
#PBS -l select=1:ppn=2

### mail alert at (b)eginning, (e)nd and (a)bortion of execution
##PBS -m bea

### send mail to the following address
##PBS -M till.schumann@rwth-aachen.de

### select the appropriate queue: short 4h, medium 32h, long 120h
##PBS -q short

### set max wallclock/CPU time
# PBS -l walltime=100:00:00
# PBS -l cput=1:00:00

### request total amount of memory
# PBS -l mem=800MB


### source openmpi
. /usr/local/mpi/openmpi/1.4.3/gcc64/bin/mpivars_openmpi-1.4.3_gcc64.sh

export LD_LIBRARY_PATH="/users/schumann/usr/szip-2.1/lib:$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH="/users/schumann/phdf5-1.8.14/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/users/schumann/usr/hdf5-1.8.14_mpi_trace/lib:$LD_LIBRARY_PATH"
export PATH="/users/schumann/Scalasca-2.1/bin:$PATH"

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo "------------------------------------------------------"
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo "------------------------------------------------------"
echo "PBS: qsub is running on $PBS_O_HOST"
echo "PBS: originating queue is $PBS_O_QUEUE"
echo "PBS: executing queue is $PBS_QUEUE"
echo "PBS: working directory is $PBS_O_WORKDIR"
echo "PBS: execution mode is $PBS_ENVIRONMENT"
echo "PBS: job identifier is $PBS_JOBID"
echo "PBS: job name is $PBS_JOBNAME"
echo "PBS: node file is $PBS_NODEFILE"
echo "PBS: current home directory is $PBS_O_HOME"
echo "PBS: PATH = $PBS_O_PATH"
echo "PBS: LD_PATH= $LD_LIBRARY_PATH"
echo "-------------------------------------------./re	-----------"

##########################################
#                                        #
#   Start job.                           #
#                                        #
##########################################

### start job from the directory it was submitted
cd $PBS_O_WORKDIR

export MPI_NODE_COUNT=1
export OMP_NUM_THREADS=2


### execute script
#CMP="scalasca -analyze -e ${PBS_JOBID}_${MPI_NODE_COUNT}_${OMP_NUM_THREADS} mpirun -np ${MPI_NODE_COUNT} ./runNESTProxy"
CMP="mpirun -np ${MPI_NODE_COUNT} ./runNESTProxy /scratch/schumann/${PBS_JOBID} nestproxyoutput_${PBS_JOBID}_${MPI_NODE_COUNT}x${OMP_NUM_THREADS}"
#CMP="mpirun -np 2 ./test_sionlib"

echo "$CMP"
printf "\n"
$CMP
