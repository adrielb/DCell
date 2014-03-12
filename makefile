DCELL_DIR=/home/abergman/projects/DCell
PETSC_DIR=/home/abergman/apps/petsc
PETSC_TMP=/home/abergman/tmp
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: run

SIM := Fibers
test: testFiberField-fiber
viz: viz-LevelSet3DView

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

sim: sim-${SIM}
run: run-${SIM}

tags:
	ctags -R  --fields=+S

cscope:
	cscope -b -q -R 

cscope-petsc:
	cd ${PETSC_DIR} && \
	cscope -b -q -R 

debug:
	export EDITOR=gvim && \
	gnome-terminal -e 'gdb FiberField/tests/fiberinit.x'

valgrind-gen-suppress:
	valgrind --leak-check=yes --gen-suppressions=all --suppressions=petscinit.supp FiberField/tests/fiberinit.x

valgrind:
	valgrind --leak-check=yes                                 \
	  --suppressions=/usr/share/openmpi/openmpi-valgrind.supp \
	  --suppressions=petscinit.supp                           \
	  ./sims/Fibers/Fibers.x

.PHONY: all build alltests rebuild opt cleanDCell sync sim run cscope cscope-petsc tags
.DEFAULT_GOAL=all
