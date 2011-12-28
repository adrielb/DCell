LOCDIR = ${DCELL_DIR}/IIM/
OBJSC = SurfaceNei.o CrossProductIIM.o SurfaceQuantities.o SurfaceDerivatives.o \
        LocalCoor.o ImmersedInterfaceContext.o ForceComponents.o Rotations.o \
        RHSUpdate.o GridToArray.o 
        #Interpolation.o 
include ${DCELL_DIR}/Common/base
LIBNAME = ImmersedInterfaceMethod

RUNOPTS = -log_summary
export TEST = corrections

all: makeLib test
	

fluidfield: ${OBJSC} makeLib 
	cd ../FluidField; make all
 
test: ${OBJSC}
	${MAKE} -C tests all

${OBJSC}: ImmersedInterfaceMethod.h IIM_private.h

#-stop_for_debugger -debugger_pause 10