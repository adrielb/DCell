#!/bin/bash

procs=8
start=0
end=0.07
JOBNAME=Fnzero2
delta=`echo "($end - $start)/($procs-1)"|bc -l`

for (( c=0; c<$procs ; c++ ))
do
    export PETSC_TMP=/data/sims/$JOBNAME/$c
    mkdir -p $PETSC_TMP
    val=`echo "$start+$delta*$c"|bc`
    echo val = $val
    
./NanoGrooves.x \
-Fa 0.1 \
-Fk 2 \
-Fn $val \
-kclip 0.1 \
-groove_width 2 \
-timax 10000 \
-CFL 0.1 \
-dtmax 1 \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-ksp_monitor \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -fieldsplit_p_ksp_max_it 4 -ksp_max_it 100 \
-log_summary -viewer_binary_skip_info \
-info $PETSC_TMP/info.log &

done
