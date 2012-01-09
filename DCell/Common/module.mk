# Common Module

LOCOBJS := LeastSq.o Serialize.o Array.o Grid.o Grid3D.o InterpolateVelocity.o \
           GlobalArrays.o DCellInit.o Memcache.o Heap.o UniqueID.o
LOCOBJS := ${addprefix ${subdirectory}/, ${LOCOBJS}}
${eval ${call make-library,${subdirectory}, ${LOCOBJS} }}

${eval ${call test-library,${subdirectory},test, Common, 2 }}
${eval ${call test-library,${subdirectory},array, Common, 1 }}
${eval ${call test-library,${subdirectory},heap, Common, 1 }}
${eval ${call test-library,${subdirectory},memcache, Common, 1 }}
${eval ${call test-library,${subdirectory},grid, Common, 1 }}
${eval ${call test-library,${subdirectory},leastsq, Common, 1 }}
${eval ${call test-library,${subdirectory},ga2da, Common, 1 }}
${eval ${call test-library,${subdirectory},uniqueid, Common, 5 }}
#${eval ${call test-library,${subdirectory},Benchmarks, Common, 1 }}

${LOCOBJS}: Common/Common.h
Common/Grid.o: Common/Grid.h
Common/Array.o: Common/Array.h
