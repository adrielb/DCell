include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: build alltests

build: ${libraries}
	@${CLINKER} -shared -o lib/${PETSC_ARCH}/libDCell.so ${OBJECTS}
	@echo Library: lib/${PETSC_ARCH}/libDCell.so

alltests: ${ALLTESTS}
	@echo Tested everything

test: build testDWorld.dworld

rebuild: cleanDCell build

opt:
	export PETSC_ARCH=gcc-opt; make rebuild

cleanDCell:
	rm -f ${OBJECTS}
	rm -f ${DCELL_DIR}/lib/${PETSC_ARCH}/lib*

sync:
	syncDCell.sh

.PHONY: all build alltests rebuild opt cleanDCell sync
