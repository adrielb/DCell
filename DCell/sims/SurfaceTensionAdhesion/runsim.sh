#!/bin/bash

export K
export fa

for (( K = 1; K <= 50; K += 1 ))
do
  for (( fa = 1; fa <= 50; fa += 1 ))
  do
    export PETSC_TMP=/scratch/RigidityAdhesion/ca.$K.$fa
    mkdir -p $PETSC_TMP
    rm -rf $PETSC_TMP/*
    qsub ./runsim.sge
#    qsub -j y -o /home/abergman/output/ca.$K.$fa -V -q all.q ./a.out -K $K -fa $fa
#    mpirun -np 1 ./a.out -K $K -fa $fa
  done
done

watchqstat.sh
