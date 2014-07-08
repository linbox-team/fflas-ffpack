#!/bin/bash

# Copyright (c) FFLAS-FFPACK
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/
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
