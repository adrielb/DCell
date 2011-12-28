#!/bin/bash

#Exit immediately if command returns non-zero exit status
set -e

export PETSC_ARCH=gcc-opt
cd $PETSC_DIR

./config/configure.py \
--with-endian=big \
--with-errorchecking=yes \
--with-debugging=no \
--with-cc=mpicc \
--with-fc=mpif90 \
--COPTFLAGS="-O3" \
--FOPTFLAGS="-O3" \
#--download-f-blas-lapack=1 


make -j 6 all
make test

#--with-mpi-dir=/share/apps/openmpi-1.4.1 \ 
#--with-mpiexec=/share/apps/openmpi-1.4.1/bin/mpiexec \

