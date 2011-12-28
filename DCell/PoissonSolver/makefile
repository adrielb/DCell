LOCDIR  = ${DCELL_DIR}/PossionSolver/
OBJSC   = IntelFFT.o #ICCG.o SpectralMethod.o 
include ${DCELL_DIR}/Common/base
LIBNAME =  PoissonSolver


#export OMP_NUM_THREADS=1

all: check


IntelFFT.o: IntelFFT.h
IntelFFT_check.o: IntelFFT.o


check: IntelFFT_check.o
	echo ${OBJSC}
	@${CLINKER} -openmp IntelFFT.o IntelFFT_check.o ${DCELL_LIB} ${GTS_LIB} ${PETSC_LIB} -lCommon
	@${MPIEXEC} -np 1 ./a.out
