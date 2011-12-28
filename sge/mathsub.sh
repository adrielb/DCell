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
# export user environment variables 
#$ -V

#set bash sesssion
# -e exit on error
# -x display command
# -u error on unset variable
set -e -x -u

#If the above is tripped, print a message before the script exits
trap 'echo "SCRIPT_FAILURE"' ERR

#limit the amount of virtual memory consumed to ~1GB
ulimit -v 1000000

export PETSC_TMP=/scratch/haptotaxis/
# USAGE
#                         $1
# qsub -t 1-10 mathsub.sh myfile.m

MFILE=$1

math < $MFILE
