#!/bin/bash
. ${DCELL_DIR}/Visualization/log2csv.sh

VIZ="python ${PWD}/${0%.sh}.py"

cd $PETSC_TMP

#cat edgelist.10000.*.csv > edgelist.10000.csv
# ls edgelist.*.csv | cut -d '.' -f 1,2 | uniq


extractParams
#column -s, -t params.csv

$VIZ

