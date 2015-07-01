#!/bin/bash

#* Copyright (c) FFLAS-FFPACK
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
CURRENT_PATH=`pwd`
HOME_PATH="${CURRENT_PATH}/.."


echo "Choose the architecture to process (`ls $HOME_PATH/Target -I CVS`):"
read ARCHI

TEST_DIR=`ls $HOME_PATH/Target/$ARCHI -I CVS -I GOTO -I ATLAS -I compilation.log | tail -n 1`
echo "Processing testing directory [$TEST_DIR]"

TEST_PATH="$HOME_PATH/Target/$ARCHI/$TEST_DIR"
PRIME=65521


if test -d ${ARCHI}; then
    echo "";
else
	mkdir ${ARCHI}
fi

f_base=fgemm
n_base=dgemm

## triangular system
f_funct=ftrsm
n_funct=dtrsm
perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${f_funct}-${PRIME}.txt > /tmp/${f_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${f_funct}-${PRIME}.txt  > /tmp/${f_funct}-GOTO.txt

perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${n_funct}-${PRIME}.txt > /tmp/${n_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${n_funct}-${PRIME}.txt  > /tmp/${n_funct}-GOTO.txt

./make_graph.sh /tmp/${f_funct}-ATLAS.txt /tmp/${f_funct}-GOTO.txt \
                  /tmp/${n_funct}-ATLAS.txt /tmp/${n_funct}-GOTO.txt \
                  "trsm" "triangular system with matrix right hand side" "${ARCHI}"

## LQUP
f_funct=lqup
n_funct=dgetrf
perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${f_funct}-${PRIME}.txt > /tmp/${f_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${f_funct}-${PRIME}.txt  > /tmp/${f_funct}-GOTO.txt

perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${n_funct}-${PRIME}.txt > /tmp/${n_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${n_funct}-${PRIME}.txt  > /tmp/${n_funct}-GOTO.txt

./make_graph.sh /tmp/${f_funct}-ATLAS.txt /tmp/${f_funct}-GOTO.txt \
                  /tmp/${n_funct}-ATLAS.txt /tmp/${n_funct}-GOTO.txt \
                  "lqup"  "matrix triangularization" "${ARCHI}"

## INVERSION
f_funct=inverse
n_funct=dgetri
perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${f_funct}-${PRIME}.txt > /tmp/${f_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${f_funct}-${PRIME}.txt  > /tmp/${f_funct}-GOTO.txt

perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${n_funct}-${PRIME}.txt > /tmp/${n_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${n_funct}-${PRIME}.txt  > /tmp/${n_funct}-GOTO.txt

./make_graph.sh /tmp/${f_funct}-ATLAS.txt /tmp/${f_funct}-GOTO.txt \
                  /tmp/${n_funct}-ATLAS.txt /tmp/${n_funct}-GOTO.txt \
                  "inversion"  "matrix inversion" "${ARCHI}"

## TRIANGULAR MATRIX INVERSION
f_funct=ftrtri
n_funct=dtrtri
perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${f_funct}-${PRIME}.txt > /tmp/${f_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${f_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${f_funct}-${PRIME}.txt  > /tmp/${f_funct}-GOTO.txt

perl make_graph_file.pl ${TEST_PATH}/ATLAS/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/ATLAS/timing-check-${n_funct}-${PRIME}.txt > /tmp/${n_funct}-ATLAS.txt
perl make_graph_file.pl ${TEST_PATH}/GOTO/timing-check-${n_base}-${PRIME}.txt \
                        ${TEST_PATH}/GOTO/timing-check-${n_funct}-${PRIME}.txt  > /tmp/${n_funct}-GOTO.txt

./make_graph.sh /tmp/${f_funct}-ATLAS.txt /tmp/${f_funct}-GOTO.txt \
                  /tmp/${n_funct}-ATLAS.txt /tmp/${n_funct}-GOTO.txt \
                  "trinversion"  "triangular matrix inversion" "${ARCHI}"
