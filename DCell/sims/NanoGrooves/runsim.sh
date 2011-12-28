#!/bin/bash

echo $HOSTNAME

export PETSC_TMP=/nas-0-0/abergman/sims
mkdir -p $PETSC_TMP

echo PETSC_TMP $PETSC_TMP

make -C ${DCELL_DIR} build
make -C ${DCELL_DIR}/sims rmTemp sim

# run the simulation
${DCELL_DIR}/sims/a.out \
-timax 1000 \
-ksp_monitor \
-ksp_rtol 1e-6 -fieldsplit_1_ksp_max_it 4 -log_summary -viewer_binary_skip_info

#-K $K -fa $fa