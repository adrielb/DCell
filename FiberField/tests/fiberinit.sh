#!/bin/bash
. ${DCELL_DIR}/Visualization/log2csv.sh

VIZ="python $PWD/FiberField/tests/fiberinit.py"

cd $PETSC_TMP

extractParams
#column -s, -t params.csv

$VIZ
