#!/bin/bash

JOBNAME=SOR

export PETSC_TMP=/data/sims/$JOBNAME/2
mkdir -p $PETSC_TMP

./NanoGrooves.x \
-Fa 1 \
-ecm 1 \
-Fk 0 \
-Fk0 30 \
-kclip 0.1 \
-groove_width 2 \
-timax 1000 \
-CFL 0.01 \
-dtmax 1 \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-ksp_monitor \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_pc_sor_its 4 \
-fieldsplit_v_fieldsplit_1_pc_sor_its 4 \
-log_summary -viewer_binary_skip_info \
-info $PETSC_TMP/info.log > $PETSC_TMP/output &

