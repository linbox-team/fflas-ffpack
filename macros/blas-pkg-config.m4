dnl  Check for BLAS
AC_DEFUN([FF_CHECK_BLAS2],
		[

		AC_MSG_CHECKING("for CBLAS (pkg-config)")

		BACKUP_CXXFLAGS=${CXXFLAGS}
		BACKUP_LIBS=${LIBS}

		BLAS_LIBS=`pkg-config --libs cblas`
		BLAS_PATH=`pkg-config --libs-only-L cblas`
		CBLAS_FLAG="-D__FFLASFFPACK_HAVE_CBLAS"

		AS_IF([test -z "${BLAS_LIBS}"],[BLAS_FOUND=false],[BLAS_FOUND=true])

		LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"
		CXXFLAGS="${BACKUP_CXXFLAGS} ${CBLAS_FLAG}"


		AC_TRY_RUN(
			[ ${CODE_CBLAS} ],
			[ blas_found="yes"
			dnl  BLAS_PATH=${BLAS_HOME}
			],
			[ blas_problem="$problem $BLAS_HOME"
			unset BLAS_LIBS ],
			[ blas_found="yes"
			blas_cross="yes"
			dnl  BLAS_PATH=${BLAS_HOME}
			])

		dnl what did we get ?
		AS_IF([ test "x$blas_found" = "xyes" ],
				[ BLAS_VENDOR="OTHER"
				AC_SUBST(BLAS_VENDOR)
				AC_SUBST(BLAS_LIBS)
				AC_SUBST(BLAS_PATH)
				AC_SUBST(CBLAS_FLAG)
				AC_DEFINE(HAVE_BLAS,1,[Define if BLAS is installed])
				AC_DEFINE(HAVE_CBLAS,1,[Define if C interface to BLAS is available])
				BLAS_FOUND=true
				AC_SUBST(BLAS_FOUND)
				AC_MSG_RESULT(found (cblas))
				],
				[AC_MSG_RESULT(not found)
				dnl  exit 1
				]
		     )


		AC_MSG_CHECKING("for CLAPACK (pkg-config)")

		LAPACK_LIBS=`pkg-config --silence-errors --libs atlas-clapack lapack`
		LAPACK_LIBS2=`pkg-config --silence-errors --libs lapack`
		AS_IF([test -z "${LAPACK_LIBS}"],[LAPACK_LIBS=${LAPACK_LIBS2}])


		AS_IF([test -z "${BLAS_PATH}"],[
		LIBS="${BACKUP_LIBS}  ${LAPACK_LIBS} ${BLAS_LIBS}" dnl llapack must be first
		],[
		LIBS="${BACKUP_LIBS} -L${BLAS_PATH} ${LAPACK_LIBS} ${BLAS_LIBS}"
		])

		dnl what did we get ?
		AC_TRY_RUN( [ ${CODE_CLAPACK} ],
				[ dgetrf_found="yes"
				clapack_found="yes"
				],
				[
				AC_TRY_RUN( [ ${CODE_CLAPACK} ],
					[ dgetrf_found="yes"],[echo "no"],[echo "cross"])
				],
				[lapack_cross="yes" ])

		AS_IF([ test "x$dgetrf_found" = "xyes" ],
				[ AC_SUBST(LAPACK_LIBS)
				AS_IF([ test "x$clapack_found" = "xyes" ],
					[AC_DEFINE(HAVE_CLAPACK,1,[Define if C interface to LAPACK is available])
					AC_MSG_RESULT( yes (clapack))
					],[
					AC_MSG_RESULT( yes (lapack))
					]
				     )
				AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is installed])
				],[
				AC_MSG_RESULT(not found (*lapack))
				])


		])
