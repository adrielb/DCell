# NanoGrooves simulation

${subdirectory}/${1}.o: ${LIBDCELL}
	@${CLINKER} ${subdirectory}/${1}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}
	@echo Simulation: sim-${1} 

RUNOPTS := \
-pls_rmin 0.1 \
-pls_rmax 0.5 \
-Fa 3 \
-Fk 1 \
-Fn 0 \
-ecm 0.033 \
-kclip 0.1 \
-groove_width 2 \
-timax 10 \
-CFL 0.1 \
-dtmax 1 \
-ksp_monitor \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -fieldsplit_p_ksp_max_it 4 -ksp_max_it 100 \
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
