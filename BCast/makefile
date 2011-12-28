PATH:=/share/apps/bin:${PATH}
LD_LIBRARY_PATH:=/share/apps/lib:${LD_LIBRARY_PATH}
CC=mpicc

all: BCast run
	
#BCast:
#	mpicc $@.c

run:
	@echo ""
	@echo "Executing........"
	mpiexec -np 4 BCast
	@echo "..............END"