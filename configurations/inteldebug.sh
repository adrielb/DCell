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

export PETSC_ARCH=inteldebug
cd $PETSC_DIR

./config/configure.py \
--with-shared=1 \
--with-vendor-compilers=intel \
--with-endian=big \
--with-blas-lapack-dir=[libmkl_lapack.so,libmkl.so,libguide.so,libpthread.so] \
--CC=icc \
--FC=ifort \
--CFLAGS="-vec-report0" \
--COPTFLAGS="-g" \
--FOPTFLAGS="-g" \
--with-errorchecking=yes \
--with-debugging=yes \
--with-mpi=0 \
--with-svn=svn \
#--LIBS="-lcheck -lglib-2.0" \
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
