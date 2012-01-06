#include ${PETSC_DIR}/bmake/common/base

myobjects = SpoolToBytes.o SpoolToBytes_check.o
CFLAGS = -vec-report0

all: run	

SpoolToBytes.o:
SpoolToBytes_check.o:


linkSpoolToTiff: SpoolToTiff.cpp
	icc -xP -O3 -ipo SpoolToTiff.cpp -o SpoolToTiff
	
runSpoolToTiff: linkSpoolToTiff
	./SpoolToTiff /home/share/Images/20070817-trouble-shoot 3
	
run: #	${CLINKER} ${myobjects}
	icc -xP -O3 SpoolToBytes.c -o SpoolToBytes
	cp SpoolToBytes /home/abergman/bin/SpoolToBytes
	
commit:
	svn ci

update:
	svn up
	
status:
	svn status
	
splint:
	splint -weak SpoolToBytes.c
	
check: SpoolToBytes_check.o
	@${CLINKER} -lcheck SpoolToBytes_check.o
	@./a.out