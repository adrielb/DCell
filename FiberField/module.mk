# FiberField module

LOCOBJS := FiberField.o SpatialBalancing.o Vertex.o Dynamics.o FluidVelocity.o Fiber.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

LOCDEPS := Common
INFOLOG:=-info ${PETSC_TMP}/info.log

${eval ${call test-library,${subdirectory},fiber,FiberField, 1 }}
${eval ${call test-library,${subdirectory},fiberinit,FiberField, 1, ${INFOLOG} -viewer_binary_skip_info}}
${eval ${call test-library,${subdirectory},ao,FiberField, 2, -log_summary }}
${eval ${call test-library,${subdirectory},balancing,FiberField, 3 }}
${eval ${call test-library,${subdirectory},mpidatatype,FiberField, 2 }}
${eval ${call test-library,${subdirectory},sameRankSendRecv,FiberField, 23 }}
${eval ${call test-library,${subdirectory},mpinei,FiberField, 9 }}
${eval ${call test-library,${subdirectory},swirl,FiberField, 9, ${INFOLOG} }}
${eval ${call test-library,${subdirectory},dynamics,FiberField, 1, ${INFOLOG} }}
${eval ${call test-library,${subdirectory},swirlfibers,FiberField, 9, ${INFOLOG} -viewer_binary_skip_info }}

${LOCOBJS}: FiberField/FiberField.h
