#!/bin/bash

cd /home/abergman/ProgramFiles/gts-0.7.6/

mkdir $PETSC_ARCH

./configure \
--prefix=/home/abergman/ProgramFiles/gts-0.7.6/$PETSC_ARCH \
CC=icc \
CFLAGS="-O3 -xP -ip -w"

make

make check

make install

