#!/bin/bash 

VIZ="python $PWD/FiberField/tests/swirl.py"

cd $PETSC_TMP

testawk() { 
  awk '
  /main/ { print $0 }
  ' info.log.* > test.csv
}

findError() {
  awk '
  /ERROR/ { print $0 }
  ' info.log.* 
}

extractX() { 
  awk '
  BEGIN { 
    OFS="," 
    print "ti,vID,x,y,z"
  }
  / time = /     { ti = $(NF) }
  / vID = /      { vID = $(NF) }
  / main.* X = / { print ti, vID, $NF }
  ' info.log.* > x.csv
}

extractParams() {
  awk '
  BEGIN { 
    OFS="," 
    print "rank,param,x,y,z"
  }
  / globalBounds.min = / { gsub(/\[|\]|,/,""); print $1,$3,$5,$6,$7 }
  / globalBounds.max = / { gsub(/\[|\]|,/,""); print $1,$3,$5,$6,$7 }
  / localBounds.min = /  { gsub(/\[|\]|,/,""); print $1,$3,$5,$6,$7 }
  / localBounds.max = /  { gsub(/\[|\]|,/,""); print $1,$3,$5,$6,$7; nextfile }
  ' info.log.* > params.csv
}

#testaw
#cat test.csv

#extractParams | tee $PETSC_TMP/params.csv
extractParams
#column -s, -t params.csv

extractX
#column -s, -t x.csv

findError

$VIZ
