#!/bin/bash

# SCALAPACK needs BLACS
# MUMPS needs SCALAPACK
# ZOLTAN needs PARMETIS
# SUPERLU_DIST needs FORTRAN

export PETSC_ARCH=iccopt
export PETSC_DIR=/home/abergman/ProgramFiles/petsc-2.3.3-p4
cd $PETSC_DIR

./config/configure.py \
--with-vendor-compilers=intel \
--with-endian=big \
--with-blas-lapack-dir=/opt/intel/mkl/9.1.023/ \
--CC=mpiicc \
--CFLAGS="-nocompchk -vec-report0 -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include" \
--COPTFLAGS="-O3 -xP -restrict -no-ipo -ip" \
--with-fortran=0 \
--with-gnu-compilers=0 \
--with-errorchecking=no \
--with-debugging=no \
--with-mpi-include=/opt/intel/mpi/2.0/include64 \
--with-mpi-lib=/opt/intel/mpi/2.0/lib64/libmpi.a \
--with-mpirun=mpiexec \
--with-svn=svn \
--with-is-color-value-type=short \
--LIBS="-lcheck -lglib-2.0" \
--download-umfpack=yes \
--download-sundials=yes \
--download-parmetis=yes \
--download-zoltan=yes \
#--download-hypre=yes \

#--download-mumps=yes \
#--download-blacs=yes \
#--download-scalapack=yes \
#--download-superlu=yes \
#--download-superlu_dist=yes


make all test
