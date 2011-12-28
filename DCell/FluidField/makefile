LOCDIR  = ${DCELL_DIR}/FluidField/
OBJSC = Solve.o FluidField.o MatrixAssembly.o Strain.o SemiLagrange.o MaskDomain.o
include ${DCELL_DIR}/Common/base
LIBNAME = FluidField

all: rmTemp makeLib test run

${OBJSC}: FluidField.h FluidField_private.h
build: rmTemp ${OBJSC} makeLib 

TEST = viewMat
ULIMIT=ulimit -v 500000; 
NP = 1
RUNOPTS = -ksp_monitor -ksp_converged_reason
#-dmmg_monitor_solution
#-fieldsplit_0_ksp_monitor -fieldsplit_0_ksp_converged_reason
#-fieldsplit_1_ksp_monitor -fieldsplit_1_ksp_converged_reason 
#-fieldsplit_0_fieldsplit_0_ksp_converged_reason # -fieldsplit_0_fieldsplit_0_ksp_monitor
#-ksp_monitor -ksp_converged_reason

testSL: rmTemp test/testSL.o
	@${CLINKER} ${PETSC_LIB} test/testSL.o
	@${MPIEXEC} -np 1 ./a.out

testViscoElastic: rmTemp test/testViscoElastic.o 
	@rm -rf a.out
	@${CLINKER} ${PETSC_LIB} test/testViscoElastic.o
	#@${MPIEXEC} -np 1 ./a.out

testCoupledStokes: rmTemp test/testCoupledStokes.o
	@${CLINKER} ${PETSC_LIB} test/testCoupledStokes.o
	@${MPIEXEC} -np 1 ./a.out -fieldsplit_1_ksp_monitor_short #-ksp_converged_reason -fieldsplit_0_ksp_converged_reason
#	 -fieldsplit_0_fieldsplit_0_ksp_converged_reason -fieldsplit_0_fieldsplit_1_ksp_converged_reason 
#  -fieldsplit_0_fieldsplit_0_ksp_monitor_short -ksp_converged_reason \
#   -fieldsplit_1_ksp_monitor_short -fieldsplit_0_ksp_monitor_short  \
# -fieldsplit_0_ksp_type cg -fieldsplit_0_pc_type sor
#-start_in_debugger
#-pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 2	

test: test/${TEST}.o makeLib 
	@${CLINKER} test/${TEST}.o ${DCELL_LIB} -lFluidField -lLevelSetMethod -lImmersedInterfaceMethod -lCommon ${PETSC_LIB}
#--mca btl_base_verbose 30 
#--mca btl_tcp_port_min_v4 2000
#-show-progress -I/home/abergman/Research/DCell/Common
#-ksp_converged_reason -ksp_rtol 1e-10 -fieldsplit_1_ksp_max_it 3000

sync:
	syncDCell.sh
	ssh levlabhn 'cd ${DCELL_DIR}/FluidField; make test run2'

testgmres: test/gmresdebug.o
	@${CLINKER} test/gmresdebug.o ${PETSC_LIB} ${GA_LIB} ${DCELL_LIB} -lCommon -lLevelSetMethod -lImmersedInterfaceMethod

#-viewer_binary_skip_info 
#-fieldsplit_1_ksp_monitor 
#-fieldsplit_0_fieldsplit_0_ksp_monitor

run2:
	${MPIEXEC} --mca btl_tcp_port_min_v4 2000 -hostfile ~/mpd2.hosts -np 1 ./a.out  ${RUNOPTS} -viewer_binary_skip_info

remote: sync
	ssh levlabhn "\
	  cd ${DCELL_DIR}/FluidField; \
	  make test run2"