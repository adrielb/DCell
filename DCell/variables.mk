#{{353 Miller, P. A. 1998;}}{{354 Mecklenburg,Robert 2005;}}
#source: Managing Projects with GNU Make, Robert Mecklenburg 2005

# $(call make-library, library-name, source-file-list)
define make-library
DCELLOBJECTS += ${2}
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

# ${call test-library,module,testprog, depends, NP}
define test-library
OBJECTS += ${1}/tests/${2}.o
ALLTESTS += test${1}-${2}
.PHONY: test${1}-${2}
${1}/tests/${2}.o: ${3:%=lib/${PETSC_ARCH}/lib%.a}
test${1}-${2}: ${1}/tests/${2}.o
	@echo "====================================================================="
	@echo Test target: $$@
	@${CLINKER} $$^ ${DCELL_LIB} ${addprefix -l, ${3} } ${PETSC_LIB}
	-rm -rf ${PETSC_TMP}/*
	@${MPIEXEC} -np ${4} ./a.out
endef

# ${call simulation,SimName,NP,runopts}
define simulation
OBJECTS += ${subdirectory}/${1}.o
.PHONY: sim-${1} run-${1}
${subdirectory}/${1}.o: ${LIBDCELL}
sim-${1}: ${subdirectory}/${1}.o
	@${CLINKER} $$^ ${DCELL_LIB} -lDCell ${PETSC_LIB}
	@echo Simulation: sim-${1}
run-${1}: sim-${1}
	@${MPIEXEC} -np ${2} ./a.out ${3}
endef

MODULES := $(subst /module.mk,,$(shell find . -name module.mk))
DCELL_INCLUDE := $(patsubst %,-I%,${MODULES})
CFLAGS += ${DCELL_INCLUDE} -I/share/apps/include/ -Wall -fPIC

GA_LIB = -L/share/apps/lib -lga
PETSC_LIB := ${PETSC_LIB} ${GA_LIB}
DCELL_LIB := -L${DCELL_DIR}/lib/${PETSC_ARCH} -Wl,-rpath,${DCELL_DIR}/lib/${PETSC_ARCH}
LIBDCELL := lib/${PETSC_ARCH}/libDCell.so

#Accumulate info from each module
DCELLOBJECTS := 
OBJECTS :=
libraries :=

#mpiexec -n 4 --bysocket --bind-to-socket --report-bindings
