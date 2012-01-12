dnl Copyright (c) FFLAS-FFPACK
dnl Written by Clément Pernet
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

AS_IF([test "x$enable_optimization" != "xno"],
[
AC_MSG_RESULT(yes)


BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

dnl  if test "x$HAVE_BLAS" = "xyes" ;then
AC_MSG_CHECKING([best threshold for Strassen-Winograd matrix multiplication])

AS_IF([test "x$HAVE_CBLAS" = "xtrue"],
[ echo "#define __FFLASFFPACK_HAVE_CBLAS 1" >> ../fflas-ffpack/fflas-ffpack-config.h])

CXXFLAGS="${BACKUP_CXXFLAGS} -I`pwd` -I`pwd`/fflas-ffpack ${BLAS_CFLAGS} ${CBLAS_FLAG}"
LIBS="${BACKUP_LIBS} ${BLAS_LIBS} "


AC_TRY_RUN([	//#define LinBoxSrcOnly
		#include <iostream>
		#include <fstream>
		//#define _LINBOX_LINBOX_CONFIG_H
		#define __FFLASFFPACK_CONFIGURATION
		#include "fflas-ffpack/config-blas.h"
		#include "fflas-ffpack/fflas-ffpack-config.h"
		#include "fflas-ffpack/field/modular-positive.h"
		#include "fflas-ffpack/fflas/fflas.h"
		#include <tests/timer.h>

		//using namespace LinBox;
		int main () {

		  FFPACK::Modular<double> F(17);
		  size_t n=1000, nmax=5000, prec=512, nbest=0, count=0;
		  Timer chrono;
		  double basetime, time;
		  bool bound=false;

		  double *A, *C;
		  A = new double[nmax*nmax];
		  C = new double[nmax*nmax];
		  for (size_t i=0; i<nmax*nmax;++i){
		    A[i]=2.;
	  	  }

		  std::ofstream outlog;
		  outlog.open("config.log", std::ofstream::out | std::ofstream::app);
		  outlog << std::endl
			 << "Threshold for finite field Strassen-Winograd matrix multiplication"
			 << std::endl;
		  do {
		    chrono.start();
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, A, n, 0., C, n, 0);
		    chrono.stop();
		    std::cout << std::endl
			      << "fgemm " << n << "x" << n << ": "
			      << chrono.usertime() << " s, "
                              << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			      << std::endl;
		    outlog << std::endl
			      << "fgemm " << n << "x" << n << ": "
			      << chrono.usertime() << " s, "
                              << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			      << std::endl;
		    basetime= chrono.usertime();
		    chrono.clear();
		    chrono.start();
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, A, n, 0., C, n, 1);
		    chrono.stop();
		    std::cout << "1Wino " << n << "x" << n << ": "
			      << chrono.usertime() << " s, "
			      << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			      << std::endl;
		    outlog << "1Wino " << n << "x" << n << ": "
			      << chrono.usertime() << " s, "
			      << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops"
			      << std::endl;
		    time= chrono.usertime();

		    if (basetime > time ){
		      count++;
		      if (count > 1){
	    		 nbest=n;
		         bound=true;
		         prec=prec>>1;
		         n-=prec;
		      }
		    }
		    else{
		      count=0;
		      if (bound)
			prec=prec>>1;
		      n+=prec;
		    }
		  } while ((prec > 64 ) && (n < nmax));

		  std::ofstream out("WinoThreshold");
		  out<<nbest;
		  out.close();

		  outlog << "defined __FFLASFFPACK_STRASSEN_OPTIMIZATION" << std::endl
			 << "defined __FFLASFFPACK_WINOTHRESHOLD to " << nbest << "" << std::endl;
	          outlog.close();

		  delete[] A;
		  delete[] C;

		  return 0;
		}
	],[
	AC_MSG_RESULT(done)
	WT="`cat WinoThreshold`"
	if test "$WT" != "0"; then
	 AC_DEFINE(STRASSEN_OPTIMIZATION,,[Define if optimized  threshold for Strassen-Winograd matrix multiplication is available])
	 AC_DEFINE_UNQUOTED(WINOTHRESHOLD, $WT, [optimized threshold for switching to strassen matrix multiplication])
	fi
	],[
	AC_MSG_RESULT(problem)
	strassen_opti="no"
	break
	],[
	AC_MSG_RESULT(cross compilation)
	strassen_opti="no"
	break
	])

],
[AC_MSG_RESULT(no)]
)

])
