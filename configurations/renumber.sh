#!/bin/bash

x=1
for i in *.png
do 
  counter=$(printf %04d $x)
  ln -sf $i img"$counter".png
  x=$(($x+1))
done
