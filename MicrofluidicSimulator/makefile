include ${DCELL_DIR}/Common/base
LOCDIR = ${DCELL_DIR}/StokesFlowProject/
OBJSC = MicrofluidicSimulator.o
#       PressureFree.o PressurePoisson.o CoupledStokes.o  DiffusionAdvection.o
AOUT = a.out

all: run

${OBJSC}: MicrofluidicSimulator.h

link: ${OBJSC}
	@echo "Linking MicrofluidicSimulator"
	${CLINKER} ${OBJSC} ${PETSC_LIB} -o ${AOUT}
#${DCELL_LIB}

run: link 
	@echo "Running StokesFlow"
	@${MPIEXEC} -np 1 ./${AOUT}

runPrint: link
	@echo "Running StokesFlow"
	@${MPIEXEC} -np 1 ./${AOUT} -ksp_converged_reason -ksp_monitor \
		-ksp_final_residual #-malloc_log -log_summary -ksp_max_it 200

check: link
	@echo "Checking"
	@${MPIEXEC} -np 1 ./${AOUT} -check #-mat_view_draw -draw_pause -1

.PHONY: link run check