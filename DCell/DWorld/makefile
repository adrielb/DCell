LOCDIR  = ${DCELL_DIR}/DWorld/
OBJSC   = DWorld.o DCell.o DCellsArray.o
include ${DCELL_DIR}/Common/base
LIBNAME = DWorld
NP ?= 1

all: makeLib test run

TEST = dworld

tests/${TEST}.o: ${OBJSC}
test: tests/${TEST}.o
	@${CLINKER} ${OBJSC} tests/${TEST}.o ${DCELL_LIB} -lLevelSetMethod -lFluidField -lImmersedInterfaceMethod -lCommon ${PETSC_LIB}

${OBJSC}: DWorld.h

Reaction.o: Reaction.h

sync:
	rsync --verbose --compress --recursive --times --checksum \
	--update --exclude '.svn' --exclude 'temp/*' --delete --ignore-times \
	/home/abergman/Research/DCell 128.220.138.90:temp/

submit:
	#mpiexec -machinefile ~/mpd.hosts -nolocal -perhost 1 -np 9 ./a.out
	mpirun -v -r ssh -f ~/mpd.hosts -genv I_MPI_DEVICE ssm -nolocal -perhost 1 -np 9 ./a.out #-snes_monitor_short -da_view -ts_monitor -log_summary
