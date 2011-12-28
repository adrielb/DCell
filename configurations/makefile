export PROGRAMFILES=/share/apps
export PREFIX=${HOME}/usr
WGET=wget -P /tmp
TARBZ2=tar -jxvC /tmp -f
TARGZ=tar  -zxvC /tmp -f

#OpenMPI http://www.open-mpi.org/software
OPENMPI=http://www.open-mpi.org/software/ompi/v1.4/downloads/
openmpi:
	if [ ! -e /tmp/${PETSC}.tar.gz ]; then ${WGET} http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-1.4.3.tar.gz; fi

#XVid Codec (http://www.xvid.org/Downloads.43.0.html)
xvid:
	wget -P /tmp http://downloads.xvid.org/downloads/xvidcore-1.1.3.tar.gz
	tar xzvf /tmp/xvidcore-1.1.3.tar.gz -C /tmp
	cd /tmp/xvidcore-1.1.3/build/generic && ./configure --prefix=${PREFIX} && make && make install
	@echo "Now recompile MPlayer"


#x264 Codec (http://www.videolan.org/developers/x264.html)
#yum install yasm
#Couldnt get yasm to work
x264:
	${WGET} ftp://ftp.videolan.org/pub/videolan/x264/snapshots/x264-snapshot-20080820-2245.tar.bz2
	tar -jxvf /tmp/x264-snapshot-20080820-2245.tar.bz2 -C /tmp
	cd /tmp/x264-snapshot-20080820-2245 && ./configure --prefix=${PREFIX} --disable-asm  && make && make install
	@echo "Now recompile MPlayer"


#MPlayer (http://www.mplayerhq.hu/DOCS/HTML-single/en)
mplayer:
	${WGET} http://www3.mplayerhq.hu/MPlayer/releases/MPlayer-1.0rc2.tar.bz2
	${TARBZ2} /tmp/MPlayer-1.0rc2.tar.bz2
	cd /tmp/MPlayer-1.0rc2 && \
	./configure --prefix=${PREFIX} --with-extralibdir=/home/abergman/usr/lib --with-extraincdir=/home/abergman/usr/include && \
	make && make install

#PETSc
#NEED TO FIGURE OUT PETSC_DIR??? USE PREFIX? INSTALL_DIR?
#http://www-unix.mcs.anl.gov/petsc/petsc-as/download/index.html
PETSC=petsc-3.2-p5
export PETSC_DIR=${PROGRAMFILES}/${PETSC}
petsc:
	if [ ! -e /tmp/${PETSC}.tar.gz ]; then ${WGET} ftp://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PETSC}.tar.gz; fi
	${TARGZ} /tmp/${PETSC}.tar.gz
	mv -u /tmp/${PETSC} ${PROGRAMFILES}
	cd ${PETSC_DIR} && petsc-gcc-debug.sh && petsc-gcc-opt.sh
#	cd ${PETSC_DIR} && intelopt.sh && inteldebug.sh

#Global Arrays
GLOBALARRAYS=ga-5-0-2
globalarrays:
	if [ -e ${PROGRAMFILES}/${GLOBALARRAYS} ]; then rm -rf ${PROGRAMFILES}/${GLOBALARRAYS}; fi
	${TARGZ} /tmp/${GLOBALARRAYS}.tgz
	cd /tmp/${GLOBALARRAYS} && \
		./configure --prefix=${PROGRAMFILES} && \
		make && make install && make checkprogs && make check MPIEXEC="mpirun -np 4"
#--with-blas="-L/share/apps/petsc-3.1-p7/gcc-debug/lib/libfblas.a" && \
#./configure CFLAGS=-g GA_COPT=-O0 --prefix=/tmp/ga-5-0-1/debug && make clean && make && make install
#		make TARGET=LINUX USE_MPI=yes && \
#		mpirun -np 2 ${PROGRAMFILES}/${GLOBALARRAYS}/global/testing/test.x
#		make TARGET=LINUX64 CC=icc FC=ifort COPT="-O3 -xP -vec-report0 -w" FOPT="-O3 -xP -vec-report0 -w" \
#		USE_MPI=yes \
#		MPI_INCLUDE=/share/apps/intel/include64 \
#		MPI_LIB=/share/apps/intel/impi/3.2.1.009/lib64/libmpi.so && \
#		echo "SUCCESS"

#/share/apps
#Execute as root
shareapps:
	su -c ' \
	mkdir /share && \
	mkdir /share/apps && \
	chown abergman:abergman /share/apps'
