dnl Copyright (c) 2012 FFLAS-FFPACK
dnl Written by Clément Pernet, Brice Boyer.
dnl This file was taken from LinBox linbox-opt.m4
dnl ========LICENCE========
dnl This file is part of the library FFLAS-FFPACK.
dnl
dnl FFLAS-FFPACK is free software: you can redistribute it and/or modify
dnl it under the terms of the  GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
dnl ========LICENCE========
dnl/






AC_DEFUN([FF_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization,
[AC_HELP_STRING([--disable-optimization], [ Disable run time optimization in FflasFpack code])])

dnl creating the optimise file unconditionally

echo "#ifndef __FFLASFFPACK_optimise_H" >  fflas-ffpack/fflas-ffpack-optimise.h
echo "#define __FFLASFFPACK_optimise_H" >> fflas-ffpack/fflas-ffpack-optimise.h
echo ""                                 >> fflas-ffpack/fflas-ffpack-optimise.h
dnl The optimise.h file has to be correcly written, so we close the #if !
echo "#endif // optimise.h"             >> fflas-ffpack/fflas-ffpack-optimise.h

AS_IF([test "x$enable_optimization" == "xyes"],
[
AC_MSG_RESULT(yes)


BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

echo "  *** OPTIMIZATION ***  "

AC_MSG_CHECKING([best threshold for Strassen-Winograd matrix multiplication])
AC_MSG_RESULT([see below])

CXXFLAGS_ALL="${BACKUP_CXXFLAGS} -I. -I.. -I`pwd` -I`pwd`/fflas-ffpack ${BLAS_CFLAGS} ${CBLAS_FLAG}"
LIBS="${BACKUP_LIBS} ${BLAS_LIBS} "
WINO=`cat optimiser/winograd.C`


dnl for Wino threshold for double
echo "  == Wino/BLAS threshold for Modular<double> == "
CXXFLAGS="${CXXFLAGS_ALL} -DFLTTYPE=Modular<double> -DOPTIMISATION_MODE"
AC_RUN_IFELSE([AC_LANG_SOURCE([${WINO}])],[
		dnl remove last line
		dnl  sed -i '$d' fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl  -i does not work on BSD sed
		sed  '$d' fflas-ffpack/fflas-ffpack-optimise.h > fflas-ffpack/fflas-ffpack-optimise.back.h ;
		mv fflas-ffpack/fflas-ffpack-optimise.back.h fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl append new definition
		cat WinoThreshold >> fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl close the file
		echo "#endif // optimise.h"  >> fflas-ffpack/fflas-ffpack-optimise.h
		dnl  echo done : `cat WinoThreshold`
		WINOT=`cat WinoThreshold |  awk  'NR==2' | awk '{print $ 3}'`
		dnl cleaning service !
		rm WinoThreshold ;
		AC_MSG_RESULT(done (${WINOT}))
		],[
		AC_MSG_RESULT(problem)
		break
		],[
		AC_MSG_RESULT(cross compilation)
		break
		])

dnl for WinoThreshold for float
echo "  == Wino/BLAS threshold for Modular<float> == "
CXXFLAGS="${CXXFLAGS_ALL} -DFLTTYPE=Modular<float> -DOPTIMISATION_MODE"
AC_RUN_IFELSE([AC_LANG_SOURCE([${WINO}])],[
		dnl remove last line
		dnl  sed -i '$ d' fflas-ffpack/fflas-ffpack-optimise.h ;
		sed  '$d' fflas-ffpack/fflas-ffpack-optimise.h > fflas-ffpack/fflas-ffpack-optimise.back.h ;
		mv fflas-ffpack/fflas-ffpack-optimise.back.h fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl append new definition
		cat WinoThreshold >> fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl close the file
		echo "#endif // optimise.h"  >> fflas-ffpack/fflas-ffpack-optimise.h
		dnl  echo done : `cat WinoThreshold`
		WINOT=`cat WinoThreshold |  awk  'NR==2' | awk '{print $ 3}'`
		dnl cleaning service !
		rm WinoThreshold ;
		AC_MSG_RESULT(done (${WINOT}))
		],[
		AC_MSG_RESULT(problem)
		break
		],[
		AC_MSG_RESULT(cross compilation)
		break
		])

dnl for Wino threshold for double
echo "  == Wino/BLAS threshold for ModularBalanced<double> == "
CXXFLAGS="${CXXFLAGS_ALL} -DFLTTYPE=ModularBalanced<double> -DOPTIMISATION_MODE"
AC_RUN_IFELSE([AC_LANG_SOURCE([${WINO}])],[
		dnl remove last line
		dnl  sed -i '$d' fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl  -i does not work on BSD sed
		sed  '$d' fflas-ffpack/fflas-ffpack-optimise.h > fflas-ffpack/fflas-ffpack-optimise.back.h ;
		mv fflas-ffpack/fflas-ffpack-optimise.back.h fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl append new definition
		cat WinoThreshold >> fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl close the file
		echo "#endif // optimise.h"  >> fflas-ffpack/fflas-ffpack-optimise.h
		dnl cleaning service !
		WINOT=`cat WinoThreshold |  awk  'NR==2' | awk '{print $ 3}'`
		dnl  echo done : `cat WinoThreshold`
		rm WinoThreshold ;
		AC_MSG_RESULT(done (${WINOT}))
		],[
		AC_MSG_RESULT(problem)
		break
		],[
		AC_MSG_RESULT(cross compilation)
		break
		])

dnl for WinoThreshold for float
echo "  == Wino/BLAS threshold for ModularBalanced<float> == "
CXXFLAGS="${CXXFLAGS_ALL} -DFLTTYPE=ModularBalanced<float> -DOPTIMISATION_MODE"
AC_RUN_IFELSE([AC_LANG_SOURCE([${WINO}])],[
		dnl remove last line
		dnl  sed -i '$ d' fflas-ffpack/fflas-ffpack-optimise.h ;
		sed  '$d' fflas-ffpack/fflas-ffpack-optimise.h > fflas-ffpack/fflas-ffpack-optimise.back.h ;
		mv fflas-ffpack/fflas-ffpack-optimise.back.h fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl append new definition
		cat WinoThreshold >> fflas-ffpack/fflas-ffpack-optimise.h ;
		dnl close the file
		echo "#endif // optimise.h"  >> fflas-ffpack/fflas-ffpack-optimise.h
		dnl  echo done : `cat WinoThreshold`
		WINOT=`cat WinoThreshold |  awk  'NR==2' | awk '{print $ 3}'`
		dnl cleaning service !
		rm WinoThreshold ;
		AC_MSG_RESULT(done (${WINOT}))
		],[
		AC_MSG_RESULT(problem)
		break
		],[
		AC_MSG_RESULT(cross compilation)
		break
		])


],
[AC_MSG_RESULT(no optimization)]
)

])
