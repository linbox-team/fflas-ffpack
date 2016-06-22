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
  FFLAS_SSEFLAG="--disable-sse"
fi

JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local
# Add path to compilers (if needed)
export PATH=$PATH:"$LOCAL_DIR/$CXX/bin"
echo $PATH
# Add specific locations (if needed)
LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$LOCAL_DIR/$CXX/lib"

# Where are blas installed (<blas_home>/lib/<blas_name>.so)
# And their name (libtotoblas)
BLAS_HOME="$LOCAL_DIR/$CXX"
BLAS_NAME=openblas

# Change these if necessary

BLAS_LIBS="-L$BLAS_HOME/lib/ -l$BLAS_NAME"
BLAS_CFLAGS=-I"$BLAS_HOME"/include

# Where is Givaro installed (using compiler CXX)
# Keep default if you did not modified PREFIX_INSTALL
GIVARO_PATH="$LOCAL_DIR/$CXX"

# Where to install fflas-ffpack binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX/$SSE"

# /!\ Warning /!\ This could be an issue if you changed
# the local installation directory
rm -rf "$PREFIX_INSTALL"/bin/fflas-ffpack* "$PREFIX_INSTALL"/include/fflas-ffpack*

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$PREFIX_INSTALL"/lib:"$GIVARO_PATH"/lib

#================#
# Setup Variables#
#================#

if [ "$CXX" == "icpc" ]; then
     distribution=`uname -m`
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
fi

#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix=$PREFIX_INSTALL --with-givaro=$GIVARO_PATH --with-blas-libs=$BLAS_LIBS --enable-optimization --enable-precompilation $FFLAS_SSEFLAG"
./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix="$PREFIX_INSTALL" --with-givaro="$GIVARO_PATH" --with-blas-libs="$BLAS_LIBS" --enable-optimization --enable-precompilation "$FFLAS_SSEFLAG"
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make prefix=$PREFIX_INSTALL install"
make PREFIX="$PREFIX_INSTALL" install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher


