# Copyright (c) 2011 Fflas-Ffpack
# This file is part of Fflas-Ffpack
# written by BB <bboyer@imag.fr>
# see COPYING for licence

AC_PREREQ(2.61)

AC_INIT([Fflas-Ffpack], [1.4.2],[ffpack-devel@googlegroups.com],[fflas-ffpack],
		[http://www.linalg.org/projects/fflas-ffpack])
AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([1.8 gnu no-dependencies -Wall -Wno-portability])
AM_CONFIG_HEADER([config.h])
AX_PREFIX_CONFIG_H(fflas-ffpack/fflas-ffpack-config.h, __FFLAFLAS)
AC_PATH_PROG(RM, rm, $FALSE)
RM="$RM -f"


AM_MAINTAINER_MODE
AM_DEP_TRACK
AM_OUTPUT_DEPENDENCY_COMMANDS

AM_ACLOCAL_INCLUDE(macros)

# work around to fix the backward compatibility issue of automake 1.10 with 1.9 (pb with MKDIR_P)
AC_SUBST([MKDIR_P])

AC_LANG([C++])

echo "-----------------------------------------------"

AC_DEBUG
AC_PROFILE
AC_WARNINGS

echo "-----------------------------------------------"
# CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
# CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

######################################################
# Try and pass different flags according to compiler #
######################################################

# disable default -g -O2 CXXFLAGS
: ${CXXFLAGS=""}

#set CXX
AC_PROG_CXX

AC_COMPILER_NAME

AC_SUBST([DEFAULT_CFLAGS])
AC_SUBST([DEBUG_CFLAGS])
AC_SUBST([TESTS_CFLAGS])

TESTS_CFLAGS="-O0"
DEBUG_CFLAGS="-g"
DEFAULT_CFLAGS="-pipe"
WARN_CFLAGS="-Wall"

if test "x$DBG" = "xyes" ; then
	DEFAULT_CFLAGS="-O0 ${DEFAULT_CFLAGS} " #those are CXXFLAGS
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DDEBUG -DFFLAFFLAS_DEBUG"
else
	DEFAULT_CFLAGS="-O2 ${DEFAULT_CFLAGS} "
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DNDEBUG -UFFLAFFLAS_DEBUG -UFFLAFFLAS_DEBUG"
fi

if test "x$PROF" = "xyes" ; then
	DEFAULT_CFLAGS="${DEFAULT_CFLAGS} -pg"
fi

if test "x$WARN" = "xyes" -o "x$WARN" = "xfull" ; then
	if test "x${CCNAM}" = "xicc" ; then
		WARN_CFLAGS="${WARN_CFLAGS} -Wcheck"
		# DEBUG_CFLAGS="-fast"
	else
		if test "x${CCNAM}" = "xgcc"  ; then
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -ansi -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long"
			fi
		else
			echo
			echo "*******************************************************"
			echo "unsupported compiler ($CCNAM). Please file a bug."
			echo "*******************************************************"
			echo
			WARN_CFLAGS="${WARN_CFLAGS}"
		fi
	fi
fi


DEFAULT_CFLAGS="${DEFAULT_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"
TESTS_CFLAGS="${TESTS_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"


AC_HEADER_STDC
AC_PROG_LIBTOOL

# newer libtool...
LT_PREREQ([2.2])
LT_INIT



echo "-----------------------------------------------"

# Machine characteristics

AC_CHECK_SIZEOF(char, 8)
AC_CHECK_SIZEOF(short, 16)
AC_CHECK_SIZEOF(int, 32)
AC_CHECK_SIZEOF(long, 32)
AC_CHECK_SIZEOF(long long, 64)
AC_CHECK_SIZEOF(__int64, 64)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

# Create some useful data types of fixed, known lengths

# We hereby assume that a character is always one byte
LINBOX_INT8="char";

case $ac_cv_sizeof_char in
  1)
     TWO_BYTES=2;
     FOUR_BYTES=4;
     EIGHT_BYTES=8;
     ;;
  8)
     TWO_BYTES=16;
     FOUR_BYTES=32;
     EIGHT_BYTES=64;
esac

case $TWO_BYTES in
  $ac_cv_sizeof_short)
     LINBOX_INT16="short";
     ;;
  $ac_cv_sizeof_int)
     LINBOX_INT16="int";
     ;;
esac

case $FOUR_BYTES in
  $ac_cv_sizeof_short)
     LINBOX_INT32="short";
     ;;
  $ac_cv_sizeof_int)
     LINBOX_INT32="int";
     ;;
  $ac_cv_sizeof_long)
     LINBOX_INT32="long";
     ;;
esac

case $EIGHT_BYTES in
  $ac_cv_sizeof_short)
     LINBOX_INT64="short";
     ;;
  $ac_cv_sizeof_int)
     LINBOX_INT64="int";
     ;;
  $ac_cv_sizeof_long)
     LINBOX_INT64="long";
     ;;
  $ac_cv_sizeof_long_long)
     LINBOX_INT64="long long";
     ;;
  $ac_cv_sizeof___int64)
     LINBOX_INT64="__int64";
     ;;
esac

AC_DEFINE_UNQUOTED(INT8, $LINBOX_INT8, Canonical 8-bit data type)
AC_DEFINE_UNQUOTED(INT16, $LINBOX_INT16, Canonical 16-bit data type)
AC_DEFINE_UNQUOTED(INT32, $LINBOX_INT32, Canonical 32-bit data type)
AC_DEFINE_UNQUOTED(INT64, $LINBOX_INT64, Canonical 64-bit data type)

echo "-----------------------------------------------"
# Feature checks
FF_MISC
# LB_DRIVER

AC_LANG_CPLUSPLUS

echo "-----------------------------------------------"

LB_CHECK_GMP(,,[
echo '*******************************************************************************'
echo ' WARNING: GMP not found!'
echo
echo ' GMP version 3.1.1 or greater with --enable-cxx is needed for this library to have Givaro. Please'
echo ' make sure GMP is installed and specify its location with the option'
echo ' --with-gmp=<prefix> when running configure.'
echo '*******************************************************************************'
])


LB_CHECK_GIVARO(,,[
echo '*******************************************************************************'
echo ' WARNING: GIVARO not found!'
echo
echo ' GIVARO version 3.4.0 or greater (<3.5) is required for some tests in this library.'
echo ' Please make sure GIVARO is installed and specify its location with the'
echo ' option --with-givaro=<prefix> when running configure.'
echo '*******************************************************************************'
])

BLAS_FOUND=false

FF_CHECK_BLAS

FF_CHECK_GOTOBLAS

FF_CHECK_GSL

if test "$BLAS_FOUND" = "false" ; then
	FF_CHECK_CBLAS
fi

if test "$BLAS_FOUND" = "false" ; then
	FF_CHECK_OTERBLAS
fi

if test "$BLAS_FOUND" = "false" ; then
	echo ''
	echo '*******************************************************************************'
	echo ' ERROR: BLAS not found!'
	echo
	echo ' BLAS routines are required for this library to compile. Please'
	echo ' make sure BLAS are installed and specify its location with the option'
	echo ' --with-blas=<lib> when running configure (or --with-cblas for atlas,...).'
	echo '*******************************************************************************'
	exit 1
fi

FF_CHECK_LAPACK

BLAS_LIBS="${BLAS_LIBS} ${LAPACK_LIBS}"
AC_SUBST(BLAS_LIBS)

# AM_CONDITIONAL(FFLAFFLAS_HAVE_BLAS,  test "x$BLAS_FOUND" != "xfalse")


# FF_BENCH

FF_OPT

FF_DOC


# if test ! -d ./benchmarks/data  ; then
	# echo "Creating data dir in benchmark" ;
	# mkdir ./benchmarks/data  ;
# fi

CXXFLAGS="${GMP_CFLAGS} ${CXXFLAGS}"
AC_SUBST(CXXFLAGS)

echo "-----------------------------------------------"

AC_CONFIG_FILES([
Makefile
macros/Makefile
fflasffpack-config
fflas-ffpack/Makefile
fflas-ffpack/fflas/Makefile
fflas-ffpack/ffpack/Makefile
fflas-ffpack/field/Makefile
utils/Makefile
doc/Makefile
tests/Makefile
benchmark/Makefile
])
AC_OUTPUT
