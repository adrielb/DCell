LOCDIR := ${DCELL_DIR}/${subdirectory}/
LOCOBJS := LeastSq.o Serialize.o Array.o Grid.o Grid3D.o InterpolateVelocity.o \
           GlobalArrays.o DCellInit.o Memcache.o Heap.o UniqueID.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${eval ${call test-library,${subdirectory},test, -lCommon, 2 }}
${eval ${call test-library,${subdirectory},array, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},heap, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},memcache, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},grid, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},leastsq, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},ga2da, -lCommon, 1 }}
${eval ${call test-library,${subdirectory},uniqueid, -lCommon, 5 }}
#${eval ${call test-library,${subdirectory},Benchmarks, -lCommon, 5 }}

${LOCOBJS}: Common/Common.h
Common/Grid.o: Common/Grid.h
Common/Array.o: Common/Array.h
