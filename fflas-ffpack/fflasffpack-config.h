#ifndef _FFLAS_FFPACK_FFLASFFPACK_CONFIG_H
#define _FFLAS_FFPACK_FFLASFFPACK_CONFIG_H 1
 
/* fflas-ffpack/fflasffpack-config.h. Generated automatically at end of configure. */
/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef __FFLAFLAS_AC_APPLE_UNIVERSAL_BUILD */

/* Define if BLAS routines are available */
#ifndef __FFLAFLAS_BLAS_AVAILABLE 
#define __FFLAFLAS_BLAS_AVAILABLE  /**/ 
#endif

/* Define if GMP is version 3.xxx */
/* #undef __FFLAFLAS_GMP_VERSION_3 */

/* Define that architecture uses big endian storage */
/* #undef __FFLAFLAS_HAVE_BIG_ENDIAN */

/* Define if BLAS is installed */
#ifndef __FFLAFLAS_HAVE_BLAS 
#define __FFLAFLAS_HAVE_BLAS  1 
#endif

/* Define if C interface to BLAS is available */
#ifndef __FFLAFLAS_HAVE_CBLAS 
#define __FFLAFLAS_HAVE_CBLAS  1 
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef __FFLAFLAS_HAVE_DLFCN_H 
#define __FFLAFLAS_HAVE_DLFCN_H  1 
#endif

/* Define to 1 if you have the <float.h> header file. */
#ifndef __FFLAFLAS_HAVE_FLOAT_H 
#define __FFLAFLAS_HAVE_FLOAT_H  1 
#endif

/* Define if GIVARO is installed */
/* #undef __FFLAFLAS_HAVE_GIVARO */

/* Define if GMP is installed */
#ifndef __FFLAFLAS_HAVE_GMP 
#define __FFLAFLAS_HAVE_GMP  1 
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef __FFLAFLAS_HAVE_INTTYPES_H 
#define __FFLAFLAS_HAVE_INTTYPES_H  1 
#endif

/* Define if lapack is available */
#ifndef __FFLAFLAS_HAVE_LAPACK 
#define __FFLAFLAS_HAVE_LAPACK  1 
#endif

/* Define to 1 if you have the <limits.h> header file. */
#ifndef __FFLAFLAS_HAVE_LIMITS_H 
#define __FFLAFLAS_HAVE_LIMITS_H  1 
#endif

/* Define that architecture uses little endian storage */
#ifndef __FFLAFLAS_HAVE_LITTLE_ENDIAN 
#define __FFLAFLAS_HAVE_LITTLE_ENDIAN  1 
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef __FFLAFLAS_HAVE_MEMORY_H 
#define __FFLAFLAS_HAVE_MEMORY_H  1 
#endif

/* Define to 1 if you have the <stddef.h> header file. */
#ifndef __FFLAFLAS_HAVE_STDDEF_H 
#define __FFLAFLAS_HAVE_STDDEF_H  1 
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef __FFLAFLAS_HAVE_STDINT_H 
#define __FFLAFLAS_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef __FFLAFLAS_HAVE_STDLIB_H 
#define __FFLAFLAS_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef __FFLAFLAS_HAVE_STRINGS_H 
#define __FFLAFLAS_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef __FFLAFLAS_HAVE_STRING_H 
#define __FFLAFLAS_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef __FFLAFLAS_HAVE_SYS_STAT_H 
#define __FFLAFLAS_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/time.h> header file. */
#ifndef __FFLAFLAS_HAVE_SYS_TIME_H 
#define __FFLAFLAS_HAVE_SYS_TIME_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef __FFLAFLAS_HAVE_SYS_TYPES_H 
#define __FFLAFLAS_HAVE_SYS_TYPES_H  1 
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef __FFLAFLAS_HAVE_UNISTD_H 
#define __FFLAFLAS_HAVE_UNISTD_H  1 
#endif

/* Canonical 16-bit data type */
#ifndef __FFLAFLAS_INT16 
#define __FFLAFLAS_INT16  short 
#endif

/* Canonical 32-bit data type */
#ifndef __FFLAFLAS_INT32 
#define __FFLAFLAS_INT32  int 
#endif

/* Canonical 64-bit data type */
#ifndef __FFLAFLAS_INT64 
#define __FFLAFLAS_INT64  long 
#endif

/* Canonical 8-bit data type */
#ifndef __FFLAFLAS_INT8 
#define __FFLAFLAS_INT8  char 
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef __FFLAFLAS_LT_OBJDIR 
#define __FFLAFLAS_LT_OBJDIR  ".libs/" 
#endif

/* Name of package */
#ifndef __FFLAFLAS_PACKAGE 
#define __FFLAFLAS_PACKAGE  "fflas-ffpack" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef __FFLAFLAS_PACKAGE_BUGREPORT 
#define __FFLAFLAS_PACKAGE_BUGREPORT  "ffpack-devel@googlegroups.com" 
#endif

/* Define to the full name of this package. */
#ifndef __FFLAFLAS_PACKAGE_NAME 
#define __FFLAFLAS_PACKAGE_NAME  "Fflas-Ffpack" 
#endif

/* Define to the full name and version of this package. */
#ifndef __FFLAFLAS_PACKAGE_STRING 
#define __FFLAFLAS_PACKAGE_STRING  "Fflas-Ffpack 1.4.0" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef __FFLAFLAS_PACKAGE_TARNAME 
#define __FFLAFLAS_PACKAGE_TARNAME  "fflas-ffpack" 
#endif

/* Define to the home page for this package. */
#ifndef __FFLAFLAS_PACKAGE_URL 
#define __FFLAFLAS_PACKAGE_URL  "http://www.linalg.org/" 
#endif

/* Define to the version of this package. */
#ifndef __FFLAFLAS_PACKAGE_VERSION 
#define __FFLAFLAS_PACKAGE_VERSION  "1.4.0" 
#endif

/* The size of `char', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF_CHAR 
#define __FFLAFLAS_SIZEOF_CHAR  1 
#endif

/* The size of `int', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF_INT 
#define __FFLAFLAS_SIZEOF_INT  4 
#endif

/* The size of `long', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF_LONG 
#define __FFLAFLAS_SIZEOF_LONG  8 
#endif

/* The size of `long long', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF_LONG_LONG 
#define __FFLAFLAS_SIZEOF_LONG_LONG  8 
#endif

/* The size of `short', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF_SHORT 
#define __FFLAFLAS_SIZEOF_SHORT  2 
#endif

/* The size of `__int64', as computed by sizeof. */
#ifndef __FFLAFLAS_SIZEOF___INT64 
#define __FFLAFLAS_SIZEOF___INT64  0 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef __FFLAFLAS_STDC_HEADERS 
#define __FFLAFLAS_STDC_HEADERS  1 
#endif

/* Define if optimized threshold for Strassen-Winograd matrix multiplication
   is available */
#ifndef __FFLAFLAS_STRASSEN_OPTIMIZATION 
#define __FFLAFLAS_STRASSEN_OPTIMIZATION  /**/ 
#endif

/* Version number of package */
#ifndef __FFLAFLAS_VERSION 
#define __FFLAFLAS_VERSION  "1.4.0" 
#endif

/* optimized threshold for switching to strassen matrix multiplication */
#ifndef __FFLAFLAS_WINOTHRESHOLD 
#define __FFLAFLAS_WINOTHRESHOLD  428 
#endif

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif
 
/* once: _FFLAS_FFPACK_FFLASFFPACK_CONFIG_H */
#endif
