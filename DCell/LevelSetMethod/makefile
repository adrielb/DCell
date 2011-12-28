LOCDIR = ${DCELL_DIR}/LevelSetMethod/
OBJSC = Curvature.o Resize.o FastMarching.o FastMarching2D.o FastMarching3D.o \
        LevelSet.o IrregularNodes.o IrregularNodes3D.o AdvectSL.o AdvectRK2.o \
        OrthogonalProjection.o Initializations.o ParticleLS.o ParticleLS2D.o ParticleLS3D.o \
        NormalDirection.o
        #DrawContour.o  
include ${DCELL_DIR}/Common/base
LIBNAME = LevelSetMethod

all: rmTemp makeLib test run

TEST = fmm2d
RUNOPTS =
#-log_summary -info

${OBJSC}: LevelSetMethod.h LSM_private.h
 
tests/${TEST}.o: LevelSetMethod.h
test: tests/${TEST}.o ${OBJSC}
	@${CLINKER} tests/${TEST}.o ${DCELL_LIB} -lLevelSetMethod -lCommon ${PETSC_LIB}
	
tests: testfmm2d

runout:
	./a.out ${RUNOPTS} > /tmp/out