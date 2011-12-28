export PETSC_ARCH=intelopt
include ${DCELL_DIR}/Common/base
LOCDIR = /home/abergman/Research/EclipseProjects/PetscLSM/ACWE

all: ImageCurvature LevelSetFromThres Curvature

ImageCurvature: ImageCurvature.o
	@${CLINKER} ImageCurvature.o ${DCELL_LIB} -lLevelSetMethod -lCommon ${PETSC_LIB} -o ImageCurvature

Curvature: Curvature.o
	@${CLINKER} Curvature.o ${DCELL_LIB} -lLevelSetMethod -lCommon ${PETSC_LIB} -o Curvature
	./Curvature -height 101 -width 84 -i /tmp/img -o /tmp/curv
	
LevelSetFromThres: LevelSetFromThres.o
	${CLINKER} LevelSetFromThres.o ${DCELL_LIB} -lLevelSetMethod -lCommon ${PETSC_LIB} -o LevelSetFromThres 
	#./LevelSetFromThres -i /tmp/thres.Real64 -o /tmp/ls.Real64 -levelset_bandwidth 10 -width 75 -height 63 > /tmp/out

FileIO.o: FileIO.h
ActiveContoursWithoutEdges.o: ActiveContoursWithoutEdges.h
ActiveContoursWithoutEdges: ActiveContoursWithoutEdges.o FileIO.o
	@echo "Linking ActiveContoursWithoutEdges"
	@${CLINKER} ActiveContoursWithoutEdges.o FileIO.o ${DCELL_LIB} -lLevelSetMethod -lCommon ${PETSC_LIB} -o ActiveContoursWithoutEdges
	./ActiveContoursWithoutEdges -image /tmp/img.Real64 -height 451 -width 421 -init /tmp/init.Real64 -o /tmp/out.Real64
