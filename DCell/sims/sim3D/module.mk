# 3D simulation

${subdirectory}/${1}.o: ${LIBDCELL}
	@${CLINKER} ${subdirectory}/${1}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}
	@echo Simulation: sim-${1} 

RUNOPTS := \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-Fa 1 \
-Fk 1 \
-Fk0 20 \
-Fn 0 \
-ecm 0.033 \
-kclip 0.1 \
-contactThres 0.4 \
-fluid_dx 0.20 \
-fluid_lens 20,20,9 \
-timax 10000 \
-CFL 0.1 \
-dtmax 0.1 \
-ksp_monitor \
-ksp_atol 1e-12 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_pc_type icc \
-fieldsplit_v_fieldsplit_1_pc_type icc \
-fieldsplit_v_fieldsplit_2_pc_type icc \
-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_2_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_0_pc_factor_levels 10 \
-fieldsplit_v_fieldsplit_1_pc_factor_levels 10 \
-fieldsplit_v_fieldsplit_2_pc_factor_levels 10 \
-fieldsplit_v_fieldsplit_0_pc_factor_fill 28.7972 \
-fieldsplit_v_fieldsplit_1_pc_factor_fill 28.7972 \
-fieldsplit_v_fieldsplit_2_pc_factor_fill 28.7972 \
-log_summary -viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log

#ASM-ICC
#-fieldsplit_v_fieldsplit_0_pc_type asm \
#-fieldsplit_v_fieldsplit_1_pc_type asm \
#-fieldsplit_v_fieldsplit_2_pc_type asm \
#-fieldsplit_v_fieldsplit_0_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_1_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_2_ksp_max_it 1 \
#-fieldsplit_v_fieldsplit_0_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_1_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_2_pc_asm_overlap 25 \
#-fieldsplit_v_fieldsplit_0_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_2_sub_ksp_type preonly \
#-fieldsplit_v_fieldsplit_0_sub_pc_type cholesky \
#-fieldsplit_v_fieldsplit_1_sub_pc_type cholesky \
#-fieldsplit_v_fieldsplit_2_sub_pc_type cholesky \
#-fieldsplit_v_fieldsplit_0_sub_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_sub_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_2_sub_pc_factor_mat_ordering_type nd \

#ICC
#-fieldsplit_v_fieldsplit_0_pc_type icc \
#-fieldsplit_v_fieldsplit_1_pc_type icc \
#-fieldsplit_v_fieldsplit_2_pc_type icc \
#-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_2_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_0_pc_factor_levels 15 \
#-fieldsplit_v_fieldsplit_1_pc_factor_levels 15 \
#-fieldsplit_v_fieldsplit_2_pc_factor_levels 15 \
#-fieldsplit_v_fieldsplit_0_pc_factor_fill 28.7972 \
#-fieldsplit_v_fieldsplit_1_pc_factor_fill 28.7972 \
#-fieldsplit_v_fieldsplit_2_pc_factor_fill 28.7972 \

#Cholesky
#-fieldsplit_v_fieldsplit_0_pc_type cholesky \
#-fieldsplit_v_fieldsplit_1_pc_type cholesky \
#-fieldsplit_v_fieldsplit_2_pc_type cholesky \
#-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
#-fieldsplit_v_fieldsplit_2_pc_factor_mat_ordering_type nd \

#-fieldsplit_v_fieldsplit_0_ksp_type preonly \
#-fieldsplit_v_fieldsplit_1_ksp_type preonly \
#-ksp_atol 1e3
#-ksp_rtol 1e-3
#-fieldsplit_p_ksp_monitor 
#-fieldsplit_1_ksp_max_it 4
#-ksp_monitor -ksp_converged_reason -viewer_binary_skip_info -log_summary
#-log_summary -ksp_view 
#-ksp_monitor_range -ksp_monitor_singular_value
#-fieldsplit_1_ksp_monitor_draw -ksp_rtol 1e-50
#-ksp_max_it 4 


${eval ${call simulation,sim3D,1,${RUNOPTS}}}
