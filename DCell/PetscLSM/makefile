LOCDIR = /home/abergman/Research/EclipseProjects/PetscLSM/
CFLAGS = -vec-report0 -nocompchk
include ${PETSC_DIR}/bmake/common/base

export C_INCLUDE = ${INCLUDE}

#${myobjects}: chkopts
#Common MyPCShell ACWE FastMarchingMethod StokesFlowProject IIM
NP=1
CurrentProject = Common
all:
#	make -C ${CurrentProject} all

#
# Block Jacobi
#
BlockJacobiProject = BlockJacobi.o

${BlockJacobiProject}:

#
# Dual Processor MPI communicator setup
#
DualProcProject = DualProc.o

${DualProcProject}:


# 
#	Subversion commands
#

commit:
	svn ci

update:
	svn up
	
status:
	svn status
	
add:
	svn add *.c
	svn add *.h
	svn add makefile
	
#
# Template Rules
#
link:
	@echo "Linking ${CurrentProject}"
	@make -C ${CurrentProject} link
	
run:
	@echo "Running ${CurrentProject}"
	@make -C ${CurrentProject} run
	
check:
	@echo "Checking ${CurrentProject}"
	@make -C ${CurrentProject} check
	
.PHONY: run link check

#
# Submit to SGE
#
submit: 
	svn ci -m "qmon submission ${CurrentProject}"
	ssh -f levlabhn ./submit.sh ${NP} a.out

mpdboot:
	mpdboot --file=/home/abergman/mpd.hosts