#!/bin/bash

END=`cat paramsfile|wc -l`

qsub -t 1-$END runsim.sge

watchqstat.sh
