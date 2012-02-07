# LevelSetMethod Module

#LOCDIR := ${DCELL_DIR}/${subdirectory}/
LOCOBJS := Curvature.o Resize.o FastMarching.o FastMarching2D.o FastMarching3D.o \
         LevelSet.o IrregularNodes.o IrregularNodes3D.o AdvectSL.o AdvectRK2.o \
         OrthogonalProjection.o Initializations.o ParticleLS.o ParticleLS2D.o \
         ParticleLS3D.o NormalDirection.o AdvectImplicit.o
         #DrawContour.o  
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

LOCDEPS := LevelSetMethod Common
${eval ${call test-library,${subdirectory},fmm2d,${LOCDEPS}, 2 }}

${LOCOBJS}: LevelSetMethod/LevelSetMethod.h LevelSetMethod/LSM_private.h
