# FiberField module

LOCOBJS := FiberField.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

LOCDEPS := Common

${eval ${call test-library,${subdirectory},fiber,Common, 1 }}
${eval ${call test-library,${subdirectory},fiberinit,FiberField, 1 }}

${LOCOBJS}: FiberField/FiberField.h


