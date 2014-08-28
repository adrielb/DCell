#{{353 Miller, P. A. 1998;}}{{354 Mecklenburg,Robert 2005;}}
#source: Managing Projects with GNU Make, Robert Mecklenburg 2005
SHELL=/bin/bash -e -o pipefail
PYTHONPATH:=${DCELL_DIR}/Visualization:${PYTHONPATH}
VPATH=${CURDIR}

${LIBDCELL}: ${libraries}
	@${CLINKER} -shared -o ${LIBDCELL} ${DCELLOBJECTS}
	@echo Library: ${LIBDCELL}

tests/%.o : ${LIBDCELL}
sims/%.o : ${LIBDCELL}
%.x : %.o
	${CLINKER} $^ ${DCELL_LIB} -lDCell ${PETSC_LIB} -o $@ 

# $(call make-library, library-name, source-file-list)
define make-library
DCELLOBJECTS += ${2}
CLEAN += ${2}
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

# ${call test-library,module,testprog, depends, NP, runopts}
define test-library
CLEAN += ${1}/tests/${2}.o ${1}/tests/${2}.x
ALLTESTS += ${1}/tests/${2}
run-${1}/tests/${2}.x: ${PETSC_TMP}/${1}/tests/${2}/stdout.log
${PETSC_TMP}/${1}/tests/${2}/stdout.log: ${1}/tests/${2}.x
	@echo "====================================================================="
	@echo Test target: ${1}/tests/${2}
	@mkdir -p $${@D}
	-find $${@D} -type f -delete
	@${MPIEXEC} -wdir $${@D} -np ${4} ${CURDIR}/${1}/tests/${2}.x ${5} \
			 > >(tee $${@D}/stdout.log) \
			2> >(tee $${@D}/stderr.log >&2)
endef

# ${call simulation,SimName,NP,runopts}
define simulation
CLEAN += ${subdirectory}/${1}.o ${subdirectory}/${1}.x
run-${subdirectory}/${1}.x: ${PETSC_TMP}/sims/${1}/${1}.log
${PETSC_TMP}/sims/${1}/stdout.log: ${subdirectory}/${1}.x ${subdirectory}/module.mk
	@echo "====================================================================="
	@echo Sim target: $$<
	@mkdir -p $${@D}
	-find $${@D} -type f -delete
	@${MPIEXEC} -wdir $${@D} -np ${2} ${CURDIR}/${subdirectory}/${1}.x ${3} \
			 > >(tee $${@D}/stdout.log) \
			2> >(tee $${@D}/stderr.log >&2)
endef

# copy good stdout.log to repo 
cache: 
	cp ${PETSC_TMP}/${TEST}/stdout.log ${CURDIR}/${TEST}.log
	-rm ${PETSC_TMP}/${TEST}/stdout.diff

# diff cache log with current log
${PETSC_TMP}/%.diff : ${PETSC_TMP}/%.log 
	@echo "Diffing " ${subst ${PETSC_TMP}/,,${@D}}
	-@CACHED=${subst ${PETSC_TMP},${CURDIR},${@D}}.log; \
	${DIFF} $$CACHED $< > >(tee $@) 2>&1

TESTSUITE=${firstword ${subst /, ,${TEST}}}
.SECONDEXPANSION:
testsuite: $${patsubst %,$${PETSC_TMP}/%/stdout.diff,$${filter $${TESTSUITE}/%,$${ALLTESTS}}} check
	@echo "Testsuite: " ${TESTSUITE}

.SECONDEXPANSION:
testall: $${patsubst %,$${PETSC_TMP}/%/stdout.diff,$${ALLTESTS}} check
	@echo "All tests"

check:
	@echo "Checking for failed tests"
	@DIFFs=`find ${PETSC_TMP} -name stdout.diff`; \
	for d in $${DIFFs}; do                        \
		if [ -s $$d ]; then                       \
			echo $$d;                             \
			FAILEDTEST=1;                         \
		fi                                        \
	done;                                         \
	if [ -z "$$FAILEDTEST" ]; then                \
		echo "ALL PASSED";                        \
	fi

# delete all files in $PETSC_TMP but not $PETSC_TMP itself
.PHONY: rmTemp
rmTemp: 
	-find ${PETSC_TMP} -not -path ${PETSC_TMP} -delete

ifndef PETSC_TMP
  ${error PETSC_TMP not set}
endif
ifndef PETSC_DIR
  ${error PETSC_DIR not set}
endif
ifndef DCELL_DIR
  ${error DCELL_DIR not set}
endif

MODULES := $(subst /module.mk,,$(shell find . -name module.mk))
DCELL_INCLUDE := $(patsubst %,-I%,${MODULES})
CFLAGS += ${DCELL_INCLUDE} -I/home/abergman/apps/ga-5-0-2/include -Wall -fPIC -Werror
#-Wextra

GA_LIB = -L/home/abergman/apps/ga-5-0-2/lib -lga
PETSC_LIB := ${PETSC_LIB} ${GA_LIB}
DCELL_LIB := -L${DCELL_DIR}/lib/${PETSC_ARCH} -Wl,-rpath,${DCELL_DIR}/lib/${PETSC_ARCH}
LIBDCELL := lib/${PETSC_ARCH}/libDCell.so

#Accumulate info from each module
DCELLOBJECTS := 
CLEAN :=
libraries :=

#mpiexec -n 4 --bysocket --bind-to-socket --report-bindings
# -find $${@D} -type f -print0 | xargs -n1000 -0 rm

#VTK building and visualization {{{
VTK_DIR=/share/apps/vtk/VTK
VTK_INCLUDE= \
-I${VTK_DIR} \
-I${VTK_DIR}/Utilities \
-I${VTK_DIR}/VolumeRendering \
-I${VTK_DIR}/Rendering \
-I${VTK_DIR}/Hybrid \
-I${VTK_DIR}/Widgets \
-I${VTK_DIR}/Rendering/Testing/Cxx \
-I${VTK_DIR}/IO \
-I${VTK_DIR}/Imaging \
-I${VTK_DIR}/Graphics \
-I${VTK_DIR}/GenericFiltering \
-I${VTK_DIR}/Filtering \
-I${VTK_DIR}/Common \
-I${VTK_DIR}/Common/Testing/Cxx \
-I${VTK_DIR}/Utilities/DICOMParser \
-I${VTK_DIR}/Utilities/vtkfreetype/include \
-I${VTK_DIR}/Utilities/vtknetcdf \
-I${VTK_DIR}/Utilities/vtkexodus2/include
VTK_LIB=-L${VTK_DIR}/bin -Wl,-rpath,${VTK_DIR}/bin                            \
-lvtkHybrid       -lvtkmetaio            -lvtkRendering  -lvtkWidgets         \
-lmpistubs        -lvtkFiltering         -lvtkImaging    -lvtksqlite          \
-lvtkalglib       -lvtkfreetype          -lvtkInfovis    -lvtksys             \
-lvtkCharts       -lvtkftgl              -lvtkIO         -lvtktiff            \
-lvtkCommon       -lvtkGenericFiltering  -lvtkjpeg       -lvtkNetCDF          \
-lvtkDICOMParser  -lvtklibxml2           -lvtkpng        -lvtkViews           \
-lvtkexoIIc       -lvtkGraphics          -lvtkproj4      -lvtkVolumeRendering \
-lvtkzlib         -lvtkverdict

# ${call vtkviz,vizName,inputFile}
define vtkviz
CLEAN += ${subdirectory}/${1}.x
${subdirectory}/${1}.x: ${subdirectory}/${1}.cxx
	@g++ -g -Wno-deprecated ${VTK_INCLUDE} ${subdirectory}/${1}.cxx ${VTK_LIB} -o $$@
viz-${1}: ${subdirectory}/${1}.x
	${subdirectory}/${1}.x ${2}
endef
#}}}
