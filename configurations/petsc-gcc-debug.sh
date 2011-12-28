#!/bin/bash

#Exit immediately if command returns non-zero exit status
set -e

export PETSC_ARCH=gcc-debug
cd $PETSC_DIR

./config/configure.py \
--with-shared-libraries=1 \
--with-endian=big \
--with-errorchecking=yes \
--with-debugging=yes \
--with-cc=mpicc \
--with-fc=mpif90 \
#--download-f-blas-lapack=1
#--download-superlu=1
#--download-hdf5=1
#--download-umfpack=

make all test
