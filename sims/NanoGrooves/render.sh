#!/bin/bash

S=1

MFILE=/home/abergman/Research/DCell/sims/NanoGrooves/render.m
RENDERSGE=/home/abergman/Research/DCell/sims/NanoGrooves/render.sge
MOVIESGE=/home/abergman/Research/DCell/sims/NanoGrooves/movie.sge

set -o errexit
for dir in "$@"
do
  export PETSC_TMP=$dir
  cd $PETSC_TMP
  echo "Rendering: " $PETSC_TMP
  END=`ls uvp.* | tail -n 1 | cut -d . -f 2` 
  JID=`qsub -terse -t 10000-$END:$S $RENDERSGE $MFILE | sed -r "s/\.(.*)//"`
  qsub -hold_jid $JID $MOVIESGE
done
