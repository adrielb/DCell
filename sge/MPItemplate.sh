#!/bin/bash
# request Bash shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N job
# change output file name
# -o output/$JOB_NAME
# make a parallel environment
# -pe mpi 4
#$ -cwd

echo "--------"
export PETSC_TMP=$TMP
echo $TMP/machines
cat $TMP/machines
echo "--------"
. /opt/intel/mpi/2.0/bin64/mpivars.sh
. /opt/intel/cce/9.0/bin/iccvars.sh
. /opt/intel/cmkl/8.0/tools/environment/mklvars64.sh

twice=`echo 2*$NSLOTS|bc`

/opt/intel/mpi/2.0/bin64/mpirun \
-r ssh \
-f $TMP/machines \
-genv I_MPI_DEVICE ssm \
-perhost 2 \
-n $twice \
/home/abergman/Research/EclipseProjects/PetscLSM/a.out -log_summary -m 2000

#then SCP results back to precis
#$SCP $TMP $IP:$REMOTE_DIR
#or keep the data on the nodes for distributed visualization
#cp $TMP $TMP/../DATA
#data needs to be copied to parent temp directory since child directory is deleted upon completion of job

# alternative MPI environmental variables
I_MPI_PIN_PROCS=all
I_MPI_FAST_COLLECTIVES=1
I_MPI_BCAST_NUM_PROCS=18
I_MPI_ALLTOALL_NUM_PROCS=18


echo "DONE"

