# Immersed Interface Method Module

LOCOBJS := SurfaceNei.o CrossProductIIM.o SurfaceQuantities.o SurfaceDerivatives.o \
        LocalCoor.o ImmersedInterfaceContext.o ForceComponents.o Rotations.o \
        RHSUpdate.o GridToArray.o 
        #Interpolation.o 
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: IIM/ImmersedInterfaceMethod.h IIM/IIM_private.h

LOCDEPS := -lIIM -lLevelSetMethod -lCommon
${eval ${call test-library,${subdirectory},corrections,${LOCDEPS}, 2 }}
