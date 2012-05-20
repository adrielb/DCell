# Immersed Interface Method Module

LOCOBJS := SurfaceIndex.o SurfaceQuantities.o SurfaceDerivatives.o \
        LocalCoor.o ImmersedInterfaceContext.o ForceComponents.o Rotations.o \
        RHSUpdate.o CrossProductIIM.o Interpolation.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: IIM/ImmersedInterfaceMethod.h IIM/IIM_private.h

LOCDEPS := IIM LevelSetMethod Common
${eval ${call test-library,${subdirectory},corrections,${LOCDEPS}, 2 }}
${eval ${call test-library,${subdirectory},surfquant,${LOCDEPS}, 1 }}
${eval ${call test-library,${subdirectory},interpolate,${LOCDEPS}, 1 }}
