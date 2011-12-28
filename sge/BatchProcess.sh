#!/bin/bash
# request Bourne again shell as shell for job
#$ -S /bin/bash
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N BatchProcess
# change output file name
#$ -o output/$JOB_NAME.$JOB_ID.$TASK_ID
# request dual processor queue
#$ -q dual.q

if [ "${HOSTNAME%%-*}" != "compute" ]
then
  echo "Do not run on head node"
  exit 1
fi

#TODO: check if '#$ -V' will export environment variables without executing '.bashrc'
. ~/.bashrc

#display the expanded value of PS4, followed by 
#the command and its expanded arguments or 
#associated word list
set -x

#Exit immediately if a simple command (see SHELL
#GRAMMAR above) exits with a non-zero status.
set -e

#If the above is tripped, print a message 
#before the script exits
trap 'echo "SCRIPT_FAILURE"' ERR

#Treat unset variables as an error when 
#performing parameter expansion
set -u

#limit the amount of virtual memory consumed to ~1GB
ulimit -v 1000000


# USAGE
# input and output file names must have 'III' string. 
#                              $1   $2           $3               $4
# qsub -t 1-10 BatchProcess.sh exec /remote/path input.III.Real64 output.III.tif

EXEC=$1
REMOTE_DIR=$2
INPUT_NAME=${3/III/$SGE_TASK_ID}
OUTPUT_NAME=${4/III/$SGE_TASK_ID}

SCP="scp -q -i /home/abergman/PrecisID/id_rsa_precis"
IP="LevLabHN@10.16.107.22"

cd $TMP
$SCP $IP:$REMOTE_DIR/$INPUT_NAME $TMP
$EXEC $INPUT_NAME $OUTPUT_NAME
$SCP $OUTPUT_NAME $IP:$REMOTE_DIR
