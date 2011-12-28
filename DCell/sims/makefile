LOCDIR  = ${DCELL_DIR}/sims/
include ${DCELL_DIR}/Common/base

all: rmTemp sim sync run

SIM = NanoGrooves

sim: ${SIM}/${SIM}.o
	-@rm a.out
	@${CLINKER} ${SIM}/${SIM}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}

RUNOPTS = \
-pls_rmin 0.5 \
-pls_rmax 0.9 \
-Fa 0.5 \
-Fk 1 \
-Fn 0 \
-kclip 0.1 \
-groove_width 2 \
-timax 1000 \
-CFL 0.1 \
-dtmax 1 \
-ksp_monitor \
-ksp_atol 1e-2 -ksp_rtol 1e-3 -fieldsplit_1_ksp_max_it 4 -ksp_max_it 100 \
-log_summary -viewer_binary_skip_info \
-info ${PETSC_TMP}/info.log

#-ksp_atol 1e3
#-ksp_rtol 1e-3 
#-fieldsplit_1_ksp_max_it 4
#-ksp_monitor -ksp_converged_reason -viewer_binary_skip_info -log_summary
#-log_summary -ksp_view 
#-ksp_monitor_range -ksp_monitor_singular_value
#-fieldsplit_1_ksp_monitor_draw -ksp_rtol 1e-50
#-ksp_max_it 4 
 
opt: 
	rm -f a.out ${SIM}/${SIM}.o
	export PETSC_ARCH=gcc-opt; \
	make -C ${DCELL_DIR}/sims rmTemp sim run;

sync:
	syncDCell.sh
