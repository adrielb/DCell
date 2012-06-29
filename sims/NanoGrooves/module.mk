# NanoGrooves simulation

${subdirectory}/${1}.o: ${LIBDCELL}
	@${CLINKER} ${subdirectory}/${1}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}
	@echo Simulation: sim-${1} 

RUNOPTS := \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-pls_sinit 128 \
-pls_edist 1.1 \
-ls_advectthres 1 \
-groove_width 2.0 \
-Fa 0.05 \
-Fk 0 \
-Fk0 20 \
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
-ksp_monitor \
-ksp_atol 1e-12 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_pc_type cholesky \
-fieldsplit_v_fieldsplit_1_pc_type cholesky \
-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
-log_summary -viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log

#-ksp_atol 1e3
#-ksp_rtol 1e-3
#-fieldsplit_p_ksp_monitor 
#-fieldsplit_1_ksp_max_it 4
#-ksp_monitor -ksp_converged_reason -viewer_binary_skip_info -log_summary
#-log_summary -ksp_view 
#-ksp_monitor_range -ksp_monitor_singular_value
#-fieldsplit_1_ksp_monitor_draw -ksp_rtol 1e-50
#-ksp_max_it 4 

${eval ${call simulation,NanoGrooves,1,${RUNOPTS}}}
