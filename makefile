PETSC_DIR=/home/abergman/apps/petsc
PETSC_TMP=/home/abergman/tmp
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: test

SIM := NanoGrooves
test: testFiberField-fiberinit
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

cscope:
	cscope -b -q -R 

cscope-petsc:
	cd ${PETSC_DIR} && \
	cscope -b -q -R 

.PHONY: all build alltests rebuild opt cleanDCell sync sim run cscope
.DEFAULT_GOAL=all
