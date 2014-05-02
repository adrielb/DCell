DCELL_DIR=/home/abergman/projects/DCell
PETSC_DIR=/home/abergman/apps/petsc
PETSC_TMP=/home/abergman/tmp
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include variables.mk
include $(addsuffix /module.mk,$(MODULES))

all: debug

SIM := Fibers
TEST := FiberField-fiberinit
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

test: test-${TEST}
valgrind: valgrind-${TEST}
debug: debug-${TEST}
sim: sim-${SIM}
run: run-${SIM}

tags:
	ctags -R  --fields=+S

cscope:
	cscope -b -q -R 

cscope-petsc:
	cd ${PETSC_DIR} && \
	cscope -b -q -R 

valgrind-gen-suppress:
	valgrind --leak-check=yes --gen-suppressions=all --suppressions=petscinit.supp FiberField/tests/fiberinit.x


.PHONY: all build alltests rebuild opt cleanDCell sync sim run cscope cscope-petsc tags
.DEFAULT_GOAL=all
