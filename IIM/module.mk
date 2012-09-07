# Immersed Interface Method Module

LOCOBJS := SurfaceQuantities.o SurfaceDerivatives.o WriteIrregularNodes.o \
        LocalCoor.o ImmersedInterfaceContext.o ForceComponents.o Rotations.o \
        IrregularNodes.o 
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: IIM/ImmersedInterfaceMethod.h IIM/IIM_private.h

LOCDEPS := IIM LevelSetMethod Common
${eval ${call test-library,${subdirectory},corrections,${LOCDEPS}, 2 }}
${eval ${call test-library,${subdirectory},surfquant,${LOCDEPS}, 1 }}
${eval ${call test-library,${subdirectory},interpolate,${LOCDEPS}, 1 }}
${eval ${call test-library,${subdirectory},surfcoor,${LOCDEPS}, 1 }}
${eval ${call test-library,${subdirectory},iimtest,${LOCDEPS}, 1 }}
${eval ${call test-library,${subdirectory},irregnodes,${LOCDEPS}, 1 }}
