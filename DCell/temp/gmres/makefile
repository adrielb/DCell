include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

all: dag.o
	-rm a.out
	@${CLINKER} dag.o ${PETSC_LIB}
	rm dag.o
	./a.out -ksp_monitor #-info -log_summary

sync:
	syncDCell.sh
	ssh levlabhn "cd ${DCELL_DIR}/temp/gmres; make all"
