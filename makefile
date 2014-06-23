export DCELL_DIR=/home/abergman/projects/DCell
PETSC_DIR=/home/abergman/apps/petsc
PETSC_TMP=/home/abergman/tmp
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: viz-sim

SIM := Fibers
TEST := FiberField-fiberinit
#viz: viz-LevelSet3DView

${LIBDCELL}: ${libraries}
	@${CLINKER} -shared -o ${LIBDCELL} ${DCELLOBJECTS}
	@echo Library: ${LIBDCELL}

build: ${LIBDCELL}

alltests: ${ALLTESTS}
	@echo Tested everything

rebuild: cleanDCell build

opt:
	export PETSC_ARCH=gcc-opt; \
	make cleanDCell; \
	make -j16 run

cleanDCell:
	rm -f ${CLEAN}
	rm -f lib/${PETSC_ARCH}/lib*

sync:
	syncDCell.sh

test: test-${TEST}
valgrind: valgrind-${TEST}
debug: debug-${TEST}
viz: viz-${TEST}
sim: sim-${SIM}
run: run-${SIM}
viz-sim: viz-sim-${SIM}
debug-sim: debug-sim-${SIM}
valgrind-sim: valgrind-sim-${SIM}

tags:
	ctags -R  --fields=+S

cscope:
	cscope -b -q -R 

OPENMPI_DIR=/home/abergman/apps/openmpi-1.8.1
tags-mpi: 
	cd ${OPENMPI_DIR} && \
	ctags -R --fields=+S ompi/mpi/c/*.c build/include/*.h

cscope-mpi:
	find "${OPENMPI_DIR}/build/include" -name "*.[ch]" > ${OPENMPI_DIR}/cscope.files
	find "${OPENMPI_DIR}/ompi/mpi/c" \
			-path "${OPENMPI_DIR}/ompi/mpi/c/profile/*" -prune -o \
			-name "*.[ch]" -print >> ${OPENMPI_DIR}/cscope.files
	cd ${OPENMPI_DIR} && \
	cscope -b -q && \
	ctags -L cscope.files

cscope-petsc:
	find ${PETSC_DIR} \
	  -path "${PETSC_DIR}/include/mpiuni/*" -prune -o \
	  -path "${PETSC_DIR}/include/finclude/*" -prune -o \
	  -path "${PETSC_DIR}/*/examples/*" -prune -o \
	  -path "${PETSC_DIR}/externalpackages/*" -prune -o \
	  -path "${PETSC_DIR}/src/docs/*" -prune -o \
	  -name "petsclog.h" -prune -o \
		-name "*.[ch]" -print > ${PETSC_DIR}/cscope.files && \
		cd ${PETSC_DIR} && \
	cscope -b -q 
# cscope.files checked by default
# -b  -  build cross-reference only, dont launch gui
# -q  -  build fast lookup tables 'cscope.in.out' and 'cscope.po.out'

valgrind-gen-suppress:
	valgrind --leak-check=yes --gen-suppressions=all --suppressions=petscinit.supp FiberField/tests/fiberinit.x

ipython:
	cd ${PETSC_TMP} && \
			tmux new ipython


.PHONY: all build alltests rebuild opt cleanDCell sync sim run cscope cscope-petsc tags
.DEFAULT_GOAL=all
