# DWorld Module

LOCOBJS := DWorld.o DCell.o DCellsArray.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${LOCOBJS}: DWorld/DWorld.h

LOCDEPS := -lDWorld -lFluidField -lIIM -lLevelSetMethod -lCommon
${eval ${call test-library,${subdirectory},dworld,${LOCDEPS}, 2 }}
