#!/bin/bash

JOBNAME=Fk
JOBID=08

export PETSC_TMP=/data/sims/$JOBNAME/$JOBID
mkdir -p $PETSC_TMP

mpiexec -wdir ${PETSC_TMP} /home/abergman/Research/DCell/sims/NanoGrooves/NanoGrooves.x \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-pls_sinit 128 \
-ls_advectthres 1 \
-groove_width 2.0 \
-Fa 0.1 \
-Fk 0 \
-Fk0 500 \
-Fn 0 \
-ecm 1 \
-kclip 0.1 \
-cell_radius 3.0 \
-contactThres 0.05 \
-adhesionRadius 0.45 \
-fluid_dx 0.20 \
-fluid_lens 35,9 \
-timax 50000 \
-CFL 0.05 \
-dtmax 0.05 \
-dtframe 1.0 \
-ksp_monitor \
-ksp_atol 1e-12 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_pc_type cholesky \
-fieldsplit_v_fieldsplit_1_pc_type cholesky \
-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
-log_summary -viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log > $PETSC_TMP/output &

#SPOOLES
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_type lu \
#-fieldsplit_v_fieldsplit_1_pc_type lu \
#-fieldsplit_v_fieldsplit_0_pc_factor_mat_solver_package spooles \
#-fieldsplit_v_fieldsplit_1_pc_factor_mat_solver_package spooles \
#-mat_spooles_symmetryflag 0 \
#-mat_spooles_symmetryflag 0 \
#-mat_spooles_ordering ND \
#-mat_spooles_ordering ND 
#mat_spooles_ordering <BestOfNDandMS,MMD,MS,ND

#superlu-dist
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_type lu \
#-fieldsplit_v_fieldsplit_1_pc_type lu \
#-fieldsplit_v_fieldsplit_0_pc_factor_mat_solver_package superlu_dist \
#-fieldsplit_v_fieldsplit_1_pc_factor_mat_solver_package superlu_dist \
#-fieldsplit_v_fieldsplit_0_mat_superlu_dist_colperm MMD_AT_PLUS_A \
#-fieldsplit_v_fieldsplit_1_mat_superlu_dist_colperm MMD_AT_PLUS_A

#eisenstat
#-fieldsplit_v_fieldsplit_0_pc_type eisenstat \
#-fieldsplit_v_fieldsplit_1_pc_type eisenstat \
#-fieldsplit_v_fieldsplit_0_pc_eisenstat_omega 1.9 \
#-fieldsplit_v_fieldsplit_1_pc_eisenstat_omega 1.9 \

#cholesky
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_type cholesky \
#-fieldsplit_v_fieldsplit_1_pc_type cholesky \
#-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \

#10) Parallel-block cholesky
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_type asm \
#-fieldsplit_v_fieldsplit_1_pc_type asm \
#-fieldsplit_v_fieldsplit_0_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_1_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_0_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_sub_pc_type cholesky \
#-fieldsplit_v_fieldsplit_1_sub_pc_type cholesky \
#-fieldsplit_v_fieldsplit_0_sub_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_sub_pc_factor_mat_ordering_type nd \

#11) Parallel-block icc
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_type asm \
#-fieldsplit_v_fieldsplit_1_pc_type asm \
#-fieldsplit_v_fieldsplit_0_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_1_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_0_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_sub_pc_type icc \
#-fieldsplit_v_fieldsplit_1_sub_pc_type icc \
#-fieldsplit_v_fieldsplit_0_sub_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_sub_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_0_sub_pc_factor_levels 10 \
#-fieldsplit_v_fieldsplit_1_sub_pc_factor_levels 10 \
 

#ICC
#-ksp_type preonly
#-ksp_max_it 1
#-pc_type icc
#-pc_factor_mat_ordering_type nd
#-pc_factor_levels 10

#SSOR
#-pc_sor_its <its>
#-pc_sor_lits <lits>
#-pc_sor_backward
#-pc_sor_local_forward
#-pc_sor_local_backward
#-pc_sor_omega
#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_pc_sor_omega 1.9 \
#-fieldsplit_v_fieldsplit_1_pc_sor_omega 1.9 \
#-fieldsplit_v_fieldsplit_0_pc_sor_symmetric \
#-fieldsplit_v_fieldsplit_1_pc_sor_symmetric \
#-fieldsplit_v_fieldsplit_0_pc_sor_local_symmetric \
#-fieldsplit_v_fieldsplit_1_pc_sor_local_symmetric \
