#!/bin/bash

FF_ATLAS=$1
FF_GOTO=$2
NUM_ATLAS=$3
NUM_GOTO=$4

FUNCTION_NAME=$5
FUNCTION_DESCR=$6
ARCH=$7

gnuplot <<EOF


set encoding iso_8859_1
set xlabel "Matrix dimension"
set ylabel "ratio"
set title "ratio of matrix multiplication / ${FUNCTION_DESCR}"

set terminal postscript  enhanced color 18
set output "${ARCH}/graph-${ARCH}-${FUNCTION_NAME}.eps"

plot [1000:10000][] "$FF_ATLAS"  using 1:(\$3/\$2) title "FFLAS/FFPACK (ATLAS)"  with lines 1 ,\
                    "$FF_GOTO"   using 1:(\$3/\$2) title "FFLAS/FFPACK (GOTO)"   with lines 2 ,\
                    "$NUM_ATLAS" using 1:(\$3/\$2) title "BLAS/LAPACK  (ATLAS)"  with lines 3 ,\
                    "$NUM_GOTO"  using 1:(\$3/\$2) title "BLAS/LAPACK  (GOTO)"   with lines 4




EOF