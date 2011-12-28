LOCDIR = ${DCELL_DIR}/Common/
OBJSC = LeastSq.o Serialize.o Array.o Grid.o Grid3D.o InterpolateVelocity.o \
		GlobalArrays.o DCellInit.o Memcache.o Heap.o UniqueID.o
include ${DCELL_DIR}/Common/base
LIBNAME = Common

all: makeLib testarray

TEST = uniqueid

test: ${OBJSC} tests/${TEST}.o
	@${CLINKER} tests/${TEST}.o ${OBJSC} ${PETSC_LIB} ${GA_LIB}
	@${MPIEXEC} -np 1 ./a.out -log_summary
#-da_processors_x 2 -da_processors_y 3 -da_processors_z 4
#-on_error_attach_debugger -display :0.0

${OBJSC}: Common.h
Serialize.o: Common.h
Grid.o: Grid.h
Array.o: Array.h
Benchmarks.o:

link: ${OBJSC}
	@echo ${PETSC_ARCH}
	@${CLINKER} ${OBJSC} ${PETSC_KSP_LIB}  ${CHECK_LIB}

linkBenchmark: Benchmarks.o
	@${CLINKER} Benchmarks.o ${PETSC_KSP_LIB}

benchmark: linkBenchmark
	@${MPIEXEC} -np 1 ./b.out -log_summary

.PHONY: link run

testtest: tests/test.o
	@${CLINKER} tests/test.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 2 ./a.out
	
testarray: tests/array.o
	@${CLINKER} tests/array.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out

testheap: tests/heap.o
	@${CLINKER} tests/heap.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out

testmemcache: tests/memcache.o
	@${CLINKER} tests/memcache.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out
	
testgrid: tests/grid.o
	@${CLINKER} tests/grid.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out

testleastsq: tests/leastsq.o
	@${CLINKER} $^ ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out
	echo $^
	
testga2da: tests/ga2da.o
	${CLINKER} tests/ga2da.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 1 ./a.out

testuniqueid: tests/uniqueid.o
	@${CLINKER} tests/uniqueid.o ${DCELL_LIB} -lCommon ${PETSC_LIB}
	@${MPIEXEC} -np 5 ./a.out
	
tests: makeLib testtest testarray testheap testmemcache testgrid testga2da testleastsq rmTemp

