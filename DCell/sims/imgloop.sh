#!/bin/bash

DIR=/scratch

for (( i = 10000; i <= 19200; i += 10 ))
do
  convert ${DIR}/g.$i.png -thumbnail 50% ${DIR}/f.$i.png
  echo "Image :" $i
done

animate '${DIR}/f.*.png'
