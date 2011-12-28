PATH:=/share/apps/bin:${PATH}
LD_LIBRARY_PATH:=/share/apps/lib:${LD_LIBRARY_PATH}
CC=mpicc

OBJS=DWorld.o DCell.o test.o

all: clean run

%.o: %.c
	mpicc -c $<

${OBJS}:
DWorld.o: DWorld.h DCell.h
DCell.o: DCell.h

link: ${OBJS}
	mpicc ${OBJS} 

run: link
	mpiexec -np 9 ./a.out
	@echo "END"

clean:
	rm -f *.o
	rm -f a.out