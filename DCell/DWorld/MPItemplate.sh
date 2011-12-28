#!/bin/bash
# request Bash shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N GrayScott
# change output file name
#$ -o /home/abergman/output/$JOB_NAME
# make a parallel environment
#$ -pe mpi 9
# export all environment variables
#$ -V
# executable is in current working directory
#$ -cwd
# request single processor queue
#$ -q all.q

echo "=================================================================="
echo $JOB_NAME.$JOB_ID.$TASK_ID
date
echo "TMP="$TMP
export PETSC_TMP=$TMP
echo $TMP/machines
cat $TMP/machines
echo "HOSTNAME="`hostname`
echo "NSLOTS="$NSLOTS
echo "--------"

twice=`echo 2*$NSLOTS|bc`

/opt/intel/mpi/2.0/bin64/mpirun \
-r ssh \
-f $TMP/machines \
-genv I_MPI_DEVICE ssm \
-perhost 1  \
-n 9 \
./a.out -log_summary -snes_monitor_short -da_view -ts_monitor


# alternative environmental variables
I_MPI_PIN_PROCS=all
I_MPI_FAST_COLLECTIVES=1
I_MPI_BCAST_NUM_PROCS=18
I_MPI_ALLTOALL_NUM_PROCS=18


echo "DONE"

