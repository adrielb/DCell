#!/bin/bash
# request Bourne again shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N ImageAnalysis
# change output file name
#$ -o output/$JOB_NAME.$JOB_ID.$TASK_ID
# request dual processor queue
#$ -q dual.q

if [ "${HOSTNAME%%-*}" != "compute" ]
then
  echo "Do not run on head node"
  exit 1
fi
. ~/.bashrc
set -x

#                 $1   $2          $3 $4
# qsub Imaging.sh exec /remote/path b output.III.tif

EXEC=$1
REMOTE_DIR=$2
MOD=$3
OUTPUT_NAME=${4/III/$SGE_TASK_ID}
IMAGE_NAME=image.$MOD.$SGE_TASK_ID.gray

SCP="scp -q -i /home/abergman/PrecisID/id_rsa_precis"
IP=" LevLabHN@10.16.107.22"

cd $TMP                                 && \
$SCP $IP:$REMOTE_DIR/$IMAGE_NAME $TMP   && \
$EXEC $IMAGE_NAME $OUTPUT_NAME          && \
$SCP $OUTPUT_NAME $IP:$REMOTE_DIR       && \
exit 0

echo "SCRIPT_FAILURE"
