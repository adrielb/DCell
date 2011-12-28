#!/bin/bash


cd $PETSC_TMP
END=`ls uvp.* | tail -n 1 | cut -d . -f 2`

qsub -t 10000-$END /home/abergman/Research/DCell/sims/SurfaceTensionAdhesion/render.sge

