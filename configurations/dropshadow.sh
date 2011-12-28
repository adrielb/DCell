#!/bin/bash

folder=shadow

if [ -d ${folder} ] 
then
  echo "Reusing ${folder}"
else
  mkdir ${folder}
fi

for d in `ls *.png`;
do
  echo $d
  convert $d \( +clone -repage +8+8 -threshold 50% -negate -modulate 75 -negate -blur 0x3 \) -transparent white +swap -mosaic ${folder}/$d
done

MakeMovie.sh
