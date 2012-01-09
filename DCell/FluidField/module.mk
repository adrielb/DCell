# FluidField Module

LOCOBJS := Solve.o FluidField.o MatrixAssembly.o Strain.o SemiLagrange.o MaskDomain.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: FluidField/FluidField.h FluidField/FluidField_private.h

LOCDEPS := -lFluidField -lIIM -lLevelSetMethod -lCommon
${eval ${call test-library,${subdirectory},viewMat,${LOCDEPS}, 2 }}


ULIMIT=ulimit -v 500000; 

#-ksp_monitor -ksp_converged_reason
#-dmmg_monitor_solution
#-fieldsplit_0_ksp_monitor -fieldsplit_0_ksp_converged_reason
#-fieldsplit_1_ksp_monitor -fieldsplit_1_ksp_converged_reason 
#-fieldsplit_0_fieldsplit_0_ksp_converged_reason # -fieldsplit_0_fieldsplit_0_ksp_monitor
#-ksp_monitor -ksp_converged_reason

# -fieldsplit_1_ksp_monitor_short -ksp_converged_reason -fieldsplit_0_ksp_converged_reason
#	 -fieldsplit_0_fieldsplit_0_ksp_converged_reason -fieldsplit_0_fieldsplit_1_ksp_converged_reason 
#  -fieldsplit_0_fieldsplit_0_ksp_monitor_short -ksp_converged_reason \
#   -fieldsplit_1_ksp_monitor_short -fieldsplit_0_ksp_monitor_short  \
# -fieldsplit_0_ksp_type cg -fieldsplit_0_pc_type sor
#-start_in_debugger
#-pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 2	

#--mca btl_base_verbose 30 
#--mca btl_tcp_port_min_v4 2000
#-show-progress -I/home/abergman/Research/DCell/Common
#-ksp_converged_reason -ksp_rtol 1e-10 -fieldsplit_1_ksp_max_it 3000

#-viewer_binary_skip_info 
#-fieldsplit_1_ksp_monitor 
#-fieldsplit_0_fieldsplit_0_ksp_monitor

#${MPIEXEC} --mca btl_tcp_port_min_v4 2000 -hostfile ~/mpd2.hosts -np 1 ./a.out  ${RUNOPTS} -viewer_binary_skip_info

remote: sync
	ssh levlabhn "\
	  cd ${DCELL_DIR}/FluidField; \
	  make test run2"