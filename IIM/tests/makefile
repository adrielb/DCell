LOCDIR = ${DCELL_DIR}/IIM/test
include ${DCELL_DIR}/Common/base

all: ${TEST}.o rmTemp
	@${CLINKER} ${TEST}.o ${DCELL_LIB} -lImmersedInterfaceMethod -lLevelSetMethod -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out -viewer_binary_skip_info #-on_error_abort
