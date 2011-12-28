LOCDIR = /home/abergman/Research/EclipseProjects/PetscLSM/MyPCShell
CFLAGS = -vec-report0 -nocompchk
NP = 1
include ${PETSC_DIR}/bmake/common/base

all: run

MyPCShell.o:

link: MyPCShell.o
	@${CLINKER} MyPCShell.o ${PETSC_SNES_LIB} -lcheck -o a.out

run: link
	@echo "Running on" ${NP}
	@${MPIRUN} -np ${NP} ./a.out