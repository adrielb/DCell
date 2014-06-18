#!/bin/bash 

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
