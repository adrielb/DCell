#{{353 Miller, P. A. 1998;}}{{354 Mecklenburg,Robert 2005;}}
#source: Managing Projects with GNU Make, Robert Mecklenburg 2005

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
ALLTESTS += test${1}-${2}
.PHONY: test${1}-${2}
${1}/tests/${2}.o: ${3:%=lib/${PETSC_ARCH}/lib%.a}
${1}/tests/${2}.x: ${1}/tests/${2}.o
	@${CLINKER} $$^ ${DCELL_LIB} ${3:%=-l%} ${PETSC_LIB} -o $$@ 
test${1}-${2}: ${1}/tests/${2}.x
	@echo "====================================================================="
	@echo Test target: $$@
	-rm -rf ${PETSC_TMP}/*
	@${MPIEXEC} -np ${4} ./${1}/tests/${2}.x ${5}
endef

# ${call simulation,SimName,NP,runopts}
define simulation
CLEAN += ${subdirectory}/${1}.o ${subdirectory}/${1}.x
.PHONY: sim-${1} run-${1}
${subdirectory}/${1}.o: ${LIBDCELL}
${subdirectory}/${1}.x: ${subdirectory}/${1}.o
	@${CLINKER} $$^ ${DCELL_LIB} -lDCell ${PETSC_LIB} -o $$@ 
sim-${1}: ${subdirectory}/${1}.x
	@echo Simulation: sim-${1}
run-${1}: sim-${1}
	-find ${PETSC_TMP}/ -print0 | xargs -n1000 -0 rm
	@${MPIEXEC} -wdir ${PETSC_TMP} -np ${2} ${CURDIR}/${subdirectory}/${1}.x ${3}
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
CLEAN :=
libraries :=

#mpiexec -n 4 --bysocket --bind-to-socket --report-bindings

#VTK building and visualization
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