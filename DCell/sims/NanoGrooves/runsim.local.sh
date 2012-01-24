#!/bin/bash

procs=8
start=0.01
end=0.1
JOBNAME=CFL
delta=`echo "($end - $start)/($procs-1)"|bc -l`

for (( c=0; c<$procs ; c++ ))
do
    export PETSC_TMP=/data/sims/$JOBNAME/$c
    mkdir -p $PETSC_TMP
    val=`echo "$start+$delta*$c"|bc`
    echo val = $val
    
./NanoGrooves.x \
-Fa 1 \
-ecm 1 \
-Fk 0 \
-Fk0 10 \
-kclip 0.1 \
-groove_width 2 \
-timax 40000 \
-CFL $val \
-dtmax 1 \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-ksp_monitor \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -fieldsplit_p_ksp_max_it 4 -ksp_max_it 100 \
-log_summary -viewer_binary_skip_info \
-info $PETSC_TMP/info.log > $PETSC_TMP/output &

done
