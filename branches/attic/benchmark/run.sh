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


TEST_NAME=`date '+%F-%Hh%M'`
ARCH_NAME=`uname -p`
CURRENT_PATH=`pwd`

if test "$ARCH_NAME"="unknown"; then
    echo "Enter the architecture name:"
    read ARCH_NAME
fi;

TARGET_PATH="${CURRENT_PATH}/Target/${ARCH_NAME}"
TEST_PATH="${TARGET_PATH}/${TEST_NAME}"
TEST_SRC_PATH="${CURRENT_PATH}/test-src"

DOING_COMPILATION="yes"
NEW_ARCH="yes"



if test -d ${TARGET_PATH}; then
    NEW_ARCH="no"
    echo "The $ARCH_NAME architecture has been already tested."
    while [ "$answer" != "yes" -a  "$answer" != "no"  ]; do
      echo -n "Do you want to perform another test? yes/no : "
      read answer
    done
    if test  "$answer" = "no"; then
	echo "You cancelled this script ! Bye Bye ..."
	exit
    else
	answer=" "
	while [ "$answer" != "yes" -a  "$answer" != "no"  ]; do
	    echo -n "Do you want to recompile the sources? yes/no : "
	    read answer
	done
	if test  "$answer" = "no"; then
	    DOING_COMPILATION="no"
	else
	    DOING_COMPILATION="yes"
	fi
    fi    
    mkdir $-P ${TEST_PATH}
else
    mkdir -p ${TARGET_PATH}
    mkdir -p ${TEST_PATH}
fi
  
echo

if test "${DOING_COMPILATION}" = "yes"; then

    echo "Launching compilation..."
    echo " errors will be redirected to ${TARGET_PATH}/compilation.log"
    echo

    if test -f ${TARGET_PATH}/compilation.log; then
	rm ${TARGET_PATH}/compilation.log
    fi
    
    export TARGET_PATH

## Launch GOTO test compilation
    BIN_PATH="${TARGET_PATH}/GOTO"
    if test -d ${BIN_PATH}; then
	rm -rf "${BIN_PATH}/*"
    else
	mkdir -p ${BIN_PATH}
    fi
    export BIN_PATH
    cd  ${CURRENT_PATH}/src/FFLAS_FFPACK
    echo "Compiling FFLAS_FFPACK with GOTO..."
    echo "Compiling FFLAS_FFPACK with GOTO..." >> ${TARGET_PATH}/compilation.log
    make -k GOTO_LINK=true && echo "compilation done" && echo
    cd  ${CURRENT_PATH}/src/BLAS_LAPACK
    echo "Compiling BLAS_LAPACK with GOTO..."
    echo "Compiling BLAS_LAPACK with GOTO..."  >> ${TARGET_PATH}/compilation.log
    make -k GOTO_LINK=true && echo "compilation done" && echo
    
    
## Launch ATLAS test compilation
    BIN_PATH="${TARGET_PATH}/ATLAS"
    if test -d ${BIN_PATH}; then
	rm -rf "${BIN_PATH}/*"
    else
	mkdir -p ${BIN_PATH}
    fi
    export BIN_PATH
    cd  ${CURRENT_PATH}/src/FFLAS_FFPACK 
    echo "Compiling FFLAS_FFPACK with ATLAS..."
    echo "Compiling FFLAS_FFPACK with ATLAS..."  >> ${TARGET_PATH}/compilation.log
   make -k ATLAS_LINK=true && echo "compilation done" && echo
    cd  ${CURRENT_PATH}/src/BLAS_LAPACK
    echo "Compiling BLAS_LAPACK with ATLAS..."
    echo "Compiling BLAS_LAPACK with ATLAS..."  >> ${TARGET_PATH}/compilation.log
    make -k ATLAS_LINK=true && echo "compilation done" && echo
else
    echo "Skipping compilation..."
    echo
fi


## launch testing phase
echo "Launching test..."
echo
export TEST_SRC_PATH
export TEST_PATH

## Run GOTO test
BIN_PATH="${TARGET_PATH}/GOTO"
TEST_PATH="${TARGET_PATH}/${TEST_NAME}/GOTO"
mkdir -p ${TEST_PATH}
export TEST_PATH
export BIN_PATH
echo "running FFLAS_FFPACK tests with GOTO..."
${TEST_SRC_PATH}/mesure-FFLAS_FFPACK.sh
echo "running BLAS_LAPACK tests with GOTO..."
${TEST_SRC_PATH}/mesure-BLAS_LAPACK.sh


## Run ATLAS test
BIN_PATH="${TARGET_PATH}/ATLAS"
TEST_PATH="${TARGET_PATH}/${TEST_NAME}/ATLAS"
mkdir -p ${TEST_PATH}
export TEST_PATH
export BIN_PATH
echo "running FFLAS_FFPACK tests with ATLAS..."
${TEST_SRC_PATH}/mesure-FFLAS_FFPACK.sh
echo "running BLAS_LAPACK tests with ATLAS..."
${TEST_SRC_PATH}/mesure-BLAS_LAPACK.sh
