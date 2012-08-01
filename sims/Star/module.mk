# NanoGrooves simulation

${subdirectory}/${1}.o: ${LIBDCELL}
	@${CLINKER} ${subdirectory}/${1}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}
	@echo Simulation: sim-${1} 

RUNOPTS := \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-pls_sinit 128 \
-pls_edist 0.5 \
-ls_advectthres 1 \
-cell_radius 3.0 \
-fluid_dx 0.5 \
-fluid_lens 24,24 \
-timax 300 \
-CFL 0.1 \
-dtmax 10 \
-tend 10 \
-ksp_monitor \
-ksp_atol 1e-12 -ksp_rtol 1e-3 -ksp_max_it 100 \
-fieldsplit_p_ksp_max_it 4 \
-fieldsplit_v_fieldsplit_0_pc_type cholesky \
-fieldsplit_v_fieldsplit_1_pc_type cholesky \
-fieldsplit_v_fieldsplit_0_pc_factor_mat_ordering_type nd \
-fieldsplit_v_fieldsplit_1_pc_factor_mat_ordering_type nd \
-log_summary -viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log

${eval ${call simulation,Star,1,${RUNOPTS}}}
