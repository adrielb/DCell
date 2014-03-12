# FiberField module

LOCOBJS := FiberField.o Vertex.o Edge.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

LOCDEPS := Common

${eval ${call test-library,${subdirectory},fiber,FiberField, 1 }}
${eval ${call test-library,${subdirectory},fiberinit,FiberField, 1 }}

${LOCOBJS}: FiberField/FiberField.h


