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
ARCH=`pwd | awk -F/ '{print $(NF-2)}'`
CXX=`pwd | awk -F/ '{print $(NF)}'`
#SSE=`pwd | awk -F/ '{print $NF}'`

JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local
# Add path to compilers (if needed)
export PATH=$PATH:/usr/local/bin:"$LOCAL_DIR/$CXX/bin"
echo "PATH = ${PATH}"

# Where are blas installed (<blas_home>/lib/<blas_name>.so)
# And their name (libtotoblas)
BLAS_HOME="$LOCAL_DIR/$CXX"
BLAS_NAME=openblas

# Change these if necessary

if [ "$ARCH" == "linbox-osx" ]; then
    BLAS_LIBS="-framework Accelerate"
else
    BLAS_LIBS="-L$BLAS_HOME/lib/ -l$BLAS_NAME"
fi

# Where to install fflas-ffpack binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX"

# Add specific locations (if needed)
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/usr/local/lib:"$LOCAL_DIR/$CXX/lib":"$PREFIX_INSTALL"/lib
echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:"$LOCAL_DIR/$CXX/lib/pkgconfig"
echo "PKG_CONFIG_PATH = ${PKG_CONFIG_PATH}"

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

# Particular case for Fedora: g++-6 <- g++
if [[ "$ARCH" == "linbox-fedora-amd64" &&  "$CXX" == "g++-6" ]]; then
    CXX="g++"
    CC=gcc
fi
if [[ "$ARCH" == "linbox-fedora-amd64" ]]; then
    BLAS_LIBS="-L/usr/lib64/atlas -lsatlas"
fi
if [ -z "$CC" ]; then
    if [[ $CXX == g++* ]]; then
        CC=`echo $CXX | sed -re 'y/++/cc/'`
    else
        CC="clang"
    fi
fi 

# /!\ Warning /!\ This could be an issue if you changed
# the local installation directory
rm -rf "$PREFIX_INSTALL"/bin/fflas-ffpack* "$PREFIX_INSTALL"/include/fflas-ffpack*

#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CC=$CC --prefix=$PREFIX_INSTALL --with-blas-libs=$BLAS_LIBS --enable-precompilation"
./autogen.sh CXX=$CXX CC=$CC --prefix="$PREFIX_INSTALL" --with-blas-libs="$BLAS_LIBS" --enable-precompilation
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make autotune"
make autotune
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make install"
make install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher


