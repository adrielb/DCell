LOCDIR  = ${DCELL_DIR}/sims/
include ${DCELL_DIR}/Common/base

all: rmTemp sim run 

SIM = haptotaxis

sim: ${SIM}.o
	-@rm a.out
	${CLINKER} ${SIM}.o ${DCELL_LIB} -lDCell ${PETSC_LIB}

sim1: SurfaceTensionAdhesionFeedback.o
	${CLINKER} SurfaceTensionAdhesionFeedback.o ${DCELL_LIB} -lDCell ${PETSC_LIB}

RUNOPTS = -viewer_binary_skip_info

#--mca btl_base_verbose 30
#--mca btl_tcp_port_min_v4 2000
#-show-progress -I/home/abergman/Research/DCell/Common
#-ksp_converged_reason -ksp_rtol 1e-10 -ksp_max_it 3000

sync:
	rsync -a --delete --exclude '*.o' /home/abergman/Research/DCell/FluidField levlabhn:/home/abergman/Research/DCell
#	rsync -a --delete --exclude '*.o' /home/abergman/Research/DCell/Common levlabhn:/home/abergman/Research/DCell

#-viewer_binary_skip_info
#-fieldsplit_1_ksp_monitor 
#-fieldsplit_0_fieldsplit_0_ksp_monitor
