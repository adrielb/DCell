include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: test
SIM := sim3D
test: testCommon-spatialindex

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
	make -j16 build

cleanDCell:
	rm -f ${CLEAN}
	rm -f lib/${PETSC_ARCH}/lib*

sync:
	syncDCell.sh

sim: sim-${SIM}
run: run-${SIM}

.PHONY: all build alltests rebuild opt cleanDCell sync sim run
