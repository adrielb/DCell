#!/bin/bash
# request Bash shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N EchoDate
# change output file name
#$ -o output/$JOB_NAME
# make a parallel environment
#  -pe mpi 9


. /opt/intel/mpi/2.0/bin64/mpivars.sh
. /opt/intel/cce/9.0/bin/iccvars.sh
. /opt/intel/cmkl/8.0/tools/environment/mklvars64.sh

#/opt/intel/mpi/2.0/bin64/mpirun -r ssh -f $TMP/machines -genv I_MPI_DEVICE ssm -perhost 1 -n $NSLOTS benchmarks/IMB_3.0/src/IMB-MPI1 
#/opt/intel/mpi/2.0/bin64/mpirun -r ssh -f $TMP/machines -genv I_MPI_DEVICE ssm -perhost 1 -n 18 benchmarks/IMB_3.0/src/IMB-MPI1

date

sleep 15
