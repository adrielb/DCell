#!/bin/bash
# request Bourne again shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job (do this on the command line)
# -N BatchProcess
# change output file name
#$ -o output/$JOB_NAME.$JOB_ID.$TASK_ID
# request dual processor queue
#$ -q dual.q
# export all environment variables
#$ -V
#$ -t 100-219

if [ "${HOSTNAME%%-*}" != "compute" ]
then
  echo "Do not run on head node"
  exit 1
fi

#echo script commands
set -x

#                                                 $1   $2           $3               $4
# qsub -t 1-10:2 -N DictyAnalysis BatchProcess.sh exec /remote/path input.III.Real64 output.III.tif

EXEC=/home/abergman/Research/DCell/VTK/a.out
DIR=/scratch/temp/
INPUT_FORM="phi.III.Real64"
OUTPUT_FORM="phi.III.jpg"
INPUT_NAME=$DIR/${INPUT_FORM/III/$SGE_TASK_ID}
OUTPUT_NAME=$DIR/${OUTPUT_FORM/III/$SGE_TASK_ID}

SCP="scp -q -i /home/abergman/PrecisID/id_rsa_precis"
IP="LevLabHN@10.16.107.22"

$EXEC $INPUT_NAME $OUTPUT_NAME && \
exit 0

echo ""
echo "SCRIPT_FAILURE"
exit 1;
