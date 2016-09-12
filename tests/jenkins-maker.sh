#!/bin/bash
# This file is part of the FFLAS-FFPACK library.
# It is distributed under the terms of the LGPL licence version 2.1 or later 
# (see COPYING)
# Created by AB - 2014/12/03
# Modified by AC - 2016/06/20
# Modified by CP - 2016/06/22
# Some influential environment variables:
#	CXX			C++ compiler command
#	CXXFLAGS	C++ compiler flags

# Note: This script is intended to be launched
# by the Jenkins web interface whenever it needs
# to compile the project.
# It is launched from the svn:trunk root directory.
# But should be stored in /<slave_jenkins_path>/makers/

SOURCE_DIRECTORY=$( cd "$( dirname "$0" )" && pwd )

#=============================#
# Change only these variables #
#=============================#
CXX=`pwd | awk -F/ '{print $(NF-2)}'`
SSE=`pwd | awk -F/ '{print $NF}'`

# Job fflas-ffpack with SSE option flag 
# by default sse is enabled
if [ "$SSE" == "withoutSSE" ]; then
  FFLAS_SSEFLAG="--disable-simd"
fi

JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local
# Add path to compilers (if needed)
export PATH=$PATH:/usr/local/bin:"$LOCAL_DIR/$CXX/bin"
echo $PATH

# Where are blas installed (<blas_home>/lib/<blas_name>.so)
# And their name (libtotoblas)
BLAS_HOME="$LOCAL_DIR/$CXX"
BLAS_NAME=openblas

# Change these if necessary

BLAS_LIBS="-L$BLAS_HOME/lib/ -l$BLAS_NAME"
BLAS_CFLAGS=-I"$BLAS_HOME"/include

# Where to install fflas-ffpack binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX/$SSE"

# Add specific locations (if needed)
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/usr/local/lib:"$LOCAL_DIR/$CXX/lib":"$PREFIX_INSTALL"/lib
echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:"$LOCAL_DIR/$CXX/lib/pkgconfig"
echo "PKG_CONFIG_PATH = ${PKG_CONFIG_PATH}"
# /!\ Warning /!\ This could be an issue if you changed
# the local installation directory
rm -rf "$PREFIX_INSTALL"/bin/fflas-ffpack* "$PREFIX_INSTALL"/include/fflas-ffpack*

#================#
# Setup Variables#
#================#

if [ "$CXX" == "icpc" ]; then
     distribution=`uname -m`
     CC=icc
     if [ "$distribution" == "i686" ]; then 	
	source /usr/local/bin/compilervars.sh ia32
     else
	source /usr/local/bin/compilervars.sh intel64
     fi
fi

# Particular case for Fedora23: g++=g++-5.3
vm_name=`uname -n | cut -d"-" -f1`
if [[ "$vm_name" == "fedora"  &&  "$CXX" == "g++-5.3" ]]; then
   CXX="g++"
   CC=gcc
fi
if [ -z "$CC" ]; then
    if [[ $CXX == g++* ]]; then
        CC=`echo $CXX | sed -re 'y/++/cc/'`
    else
        CC="clang"
    fi
fi 
#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CC=$CC --prefix=$PREFIX_INSTALL --with-blas-libs=$BLAS_LIBS --enable-optimization --enable-precompilation $FFLAS_SSEFLAG"
./autogen.sh CXX=$CXX CC=$CC --prefix="$PREFIX_INSTALL" --with-blas-libs="$BLAS_LIBS" --enable-optimization --enable-precompilation "$FFLAS_SSEFLAG"
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make autotune"
make autotune
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make install"
make install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher


