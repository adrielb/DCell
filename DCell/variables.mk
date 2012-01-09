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

define test-library
ALLTESTS += test${1}.${2}
.PHONY: test${1}.${2}
${1}/tests/${2}.o: lib/${PETSC_ARCH}/lib${1}.a
test${1}.${2}: ${1}/tests/${2}.o
	@echo "====================================================================="
	@echo Test target: $$@
	@${CLINKER} $$^ ${DCELL_LIB} ${3} ${PETSC_LIB}
	-rm -rf ${PETSC_TMP}/*
	@${MPIEXEC} -np ${4} ./a.out
endef

#export MODULES=Common LevelSetMethod IIM FluidField DWorld
MODULES := $(subst /module.mk,,$(shell find . -name module.mk))
DCELL_INCLUDE := $(patsubst %,-I%,${MODULES})
CFLAGS += ${DCELL_INCLUDE} -I/share/apps/include/ -Wall

GA_LIB = -L/share/apps/lib -lga
PETSC_LIB := ${PETSC_LIB} ${GA_LIB}
DCELL_LIB = -L${DCELL_DIR}/lib/${PETSC_ARCH} -Wl,-rpath,${DCELL_DIR}/lib/${PETSC_ARCH}

#Accumulate info from each module
OBJECTS := 
libraries :=
