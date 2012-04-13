LOCDIR = ${DCELL_DIR}/FiberField/
OBJSC = FiberField.o Vertex.o
include ${DCELL_DIR}/Common/base
LIBNAME = FiberField

all: makeLib test run
TEST=fiber

${OBJSC}: FiberField.h
 
tests/${TEST}.o: FiberField.h
test: tests/${TEST}.o ${OBJSC}
	@${CLINKER} tests/${TEST}.o ${DCELL_LIB} -lFiberField -lCommon ${PETSC_LIB}
