#!/bin/bash

VIZ="python $PWD/FiberField/tests/dynamics.py"

cd $PETSC_TMP

$VIZ
