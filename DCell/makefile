include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

#source: Managing Projects with GNU Make, Robert Mecklenburg 2004
# $(call make-library, library-name, source-file-list)
define make-library
OBJECTS += ${2}
libname = ./lib/${PETSC_ARCH}/lib${1}.a
libraries += ./lib/${PETSC_ARCH}/lib${1}.a
./lib/${PETSC_ARCH}/lib${1}.a: ${2}
	@echo Library: $$@
	@${AR} ${AR_FLAGS} $$@ $$^
endef

# $(subdirectory)
subdirectory = $(patsubst %/module.mk,%, \
                 $(word \
      $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))


#export MODULES=Common LevelSetMethod IIM FluidField DWorld
MODULES := $(subst /module.mk,,$(shell find . -name module.mk))
DCELL_INCLUDE := $(patsubst %,-I%,${MODULES})
CFLAGS += ${DCELL_INCLUDE} -I/share/apps/include/ -Wall

GA_LIB = -L/share/apps/lib -lga -lfblas
PETSC_LIB := ${PETSC_LIB} ${GA_LIB}
DCELL_LIB = -L${DCELL_DIR}/lib/${PETSC_ARCH} -Wl,-rpath,${DCELL_DIR}/lib/${PETSC_ARCH}

#Collect info from each module
OBJECTS := 
libraries :=

include $(addsuffix /module.mk,$(MODULES))

.PHONY: all
all: ${libraries}
	@echo MYSUB = ${MYSUB}
	@echo libraries = ${libraries}
	@echo libname = ${libname}
	@echo OBJECTS = ${OBJECTS}

rebuild:
	@clear
	@make -C Common cleanLib > /dev/null
	@echo "PETSC_ARCH = " ${PETSC_ARCH}
	@START_TIME=`date +%s`; \
	for d in ${MODULES}; do \
		echo Making $$d; \
		${MAKE} -C $$d clean makeLib > /dev/null; \
		echo ; echo; \
	done; \
	END_TIME=`date +%s.%N`; \
	echo "BUILD TIME: " `echo "$$END_TIME - $$START_TIME" | bc` " sec"\

build:
	@for d in ${MODULES}; do \
		make -C $$d makeLib; \
	done;

opt:
	export PETSC_ARCH=gcc-opt; make rebuild

cleanLib:
	rm -f ${DCELL_DIR}/lib/${PETSC_ARCH}/lib*

makeLib: ${OBJSC}
#	@${CLINKER} -shared -o ${DCELL_DIR}/lib/${PETSC_ARCH}/lib${LIBNAME}.so ${OBJSC}
	@${AR} ${AR_FLAGS}     ${DCELL_DIR}/lib/${PETSC_ARCH}/lib${LIBNAME}.a  ${OBJSC}
	@${AR} ${AR_FLAGS}     ${DCELL_DIR}/lib/${PETSC_ARCH}/libDCell.a       ${OBJSC}
	@echo        "Library: ${DCELL_DIR}/lib/${PETSC_ARCH}/lib${LIBNAME}   [${OBJSC}]"

rmTemp:
	-rm -rf ${PETSC_TMP}/*

run:
	@${MPIEXEC} -np ${NP} ./a.out ${RUNOPTS}
NP ?= 1

dcellshared:
	${CLINKER} -shared -o ${DCELL_DIR}/lib/${PETSC_ARCH}/libDCell.so ${patsubst %, ${DCELL_DIR}/%/*.o, ${MODULES} }
