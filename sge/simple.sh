#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/sh
# merge stderr with stdout
#$ -j y
# assign a name to submitted job
#$ -N StreamBenchmark
# change output file name
#$ -o output/$JOB_NAME.$TASK_ID

#
# set temp directory
MYTEMPDIR=/state/partition1/temp

#
# clear out the temporary directory
rm -rf $MYTEMPDIR/*

#
# compile stream benchmark
icc -DN=$SGE_TASK_ID -fast -O ~/stream/stream.c -o $MYTEMPDIR/stream -w -vec-report0 2>/dev/null

#if [ $? ]
#then echo "ICC failed"
#fi

#
# execute stream
$MYTEMPDIR/stream
