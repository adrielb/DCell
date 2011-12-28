#!/bin/sh
# request Bourne shell as shell for job
#$ -S /bin/sh
# merge stderr with stdout
#$ -j y

rm /state/partition1/temp/*
