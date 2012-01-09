include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: sim

${LIBDCELL}: ${libraries}
	@${CLINKER} -shared -o ${LIBDCELL} ${DCELLOBJECTS}
	@echo Library: ${LIBDCELL}

build: ${LIBDCELL}

alltests: ${ALLTESTS}
	@echo Tested everything

test: build testDWorld.dworld

rebuild: cleanDCell build

opt:
	export PETSC_ARCH=gcc-opt; make rebuild

cleanDCell:
	rm -f ${OBJECTS}
	rm -f lib/${PETSC_ARCH}/lib*

sync:
	syncDCell.sh

SIM := NanoGrooves
sim: sim-${SIM}
run: run-${SIM}

.PHONY: all build alltests rebuild opt cleanDCell sync sim run
