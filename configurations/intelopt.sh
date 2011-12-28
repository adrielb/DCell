#!/bin/bash

#Exit immediately if command returns non-zero exit status
set -e

# need to undefine stricmp in ./bmake/intelopt/petscconfig.h !!!
# bin64/mpiexec defined for MPIEXEC, perhaps specify explicitly???
# specify shared libraries for blas-lapack-pdepoisson
# using "-align" optimization flag? 
# to get IPO: ln -s /usr/lib/gcc/.../libstdc++_shared.so /opt/intel/cce/.../lib

# SCALAPACK needs BLACS
# MUMPS needs SCALAPACK
# ZOLTAN needs PARMETIS
# SUPERLU_DIST needs FORTRAN

export PETSC_ARCH=intelopt
cd $PETSC_DIR

./config/configure.py \
--with-shared=1 \
--with-vendor-compilers=intel \
--with-endian=big \
--with-blas-lapack-lib=[libmkl_lapack.so,libmkl.so,libguide.so,libpthread.so] \
--CC=icc \
--FC=ifort \
--CFLAGS="-vec-report0" \
--COPTFLAGS="-O3 -xP -restrict" \
--FOPTFLAGS="-O3 -xP" \
--with-errorchecking=yes \
--with-debugging=no \
--with-mpi-shared=1 \
--with-mpi-include=$I_MPI_ROOT/include64 \
--with-mpi-lib=[$I_MPI_ROOT/lib64/libmpi.so,$I_MPI_ROOT/lib64/libmpiif.so] \
--with-mpirun=$I_MPI_ROOT/bin64/mpiexec \
--with-svn=svn \
#--with-is-color-value-type=short \
#--with-sundials=1 \
#--with-sundials-dir=$PETSC_DIR/externalpackages/sundials-2.3.0/intelopt \
#--with-umfpack=1 \
#--with-umfpack-dir=$PETSC_DIR/externalpackages/UMFPACKv4.3/intelopt \
#--download-umfpack=yes \
#--download-sundials=yes \
#--download-blacs=yes \
#--download-scalapack=yes \
#--download-superlu=yes \
#--download-spai=yes \
#--download-zoltan=yes \
#--download-mumps=yes \
#--download-parmetis=yes \
#--download-superlu_dist=yes \
#--download-hypre=yes \
#--with-gnu-compilers=0 \

make all test

#SUNDIALS config
#./configure --prefix=$PETSC_DIR/externalpackages/sundials-2.3.0/intelopt CC="icc" --with-cflags="-vec-report0 -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -fPIC -O3 -xP -restrict " --without-mpicc --with-mpi-incdir="/opt/intel/mpi/2.0/include64" --with-mpi-libdir="/opt/intel/mpi/2.0/lib64" --with-mpi-libs="-lmpigf -lmpiif -lmpi.shm -lmpi -lmpi.ssm -ltvmpi -lmpigi" --without-mpif77 --enable-examples --disable-cvodes --disable-ida --disable-kinsol --disable-f77 --disable-libtool-lock --enable-shared

#--with-mpi-lib=[/opt/intel/mpi/2.0/lib64/libmpigf.so,/opt/intel/mpi/2.0/lib64/libmpiif.so,/opt/intel/mpi/2.0/lib64/libmpi.shm.so,/opt/intel/mpi/2.0/lib64/libmpi.so,/opt/intel/mpi/2.0/lib64/libmpi.ssm.so,/opt/intel/mpi/2.0/lib64/libtvmpi.so,/opt/intel/mpi/2.0/lib64/libmpigi.a] 
