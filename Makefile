VERSION=1.3.1

FFLAS_FFPACK_DIR=/home/pernet/Logiciels/fflas-ffpack
LINBOX_DIR=/home/pernet/Logiciels/linbox


dist:
	mkdir /tmp/fflas-ffpack-${VERSION}
	mkdir /tmp/fflas-ffpack-${VERSION}/tests /tmp/fflas-ffpack-${VERSION}/include /tmp/fflas-ffpack-${VERSION}/include/fflas-ffpack /tmp/fflas-ffpack-${VERSION}/benchmark /tmp/fflas-ffpack-${VERSION}/benchmark/graph /tmp/fflas-ffpack-${VERSION}/benchmark/html /tmp/fflas-ffpack-${VERSION}/benchmark/src /tmp/fflas-ffpack-${VERSION}/benchmark/test-src /tmp/fflas-ffpack-${VERSION}/benchmark/src/BLAS_LAPACK /tmp/fflas-ffpack-${VERSION}/benchmark/src/FFLAS_FFPACK /tmp/fflas-ffpack-${VERSION}/benchmark/src/BLOCKING
	cp ${FFLAS_FFPACK_DIR}/{AUTHORS,ChangeLog,COPYING,INSTALL,README} /tmp/fflas-ffpack-${VERSION}/
	cp ${FFLAS_FFPACK_DIR}/tests/{timer.h,timer.C,Matio.h,test-fgemm.C,test-fgemv.C,test-lqup.C,test-charpoly.C,test-compressQ.C,test-frobenius.C,test-fsquare.C,test-det.C,test-invert.C,test-krylov-elim.C,test-rank.C,test-ftrsm.C,testeur_fgemm.C,testeur_lqup.C,testeur_ftrsm.C,dense_generator.C}  /tmp/fflas-ffpack-${VERSION}/tests/
	cp ${FFLAS_FFPACK_DIR}/tests/Makefile.template /tmp/fflas-ffpack-${VERSION}/tests/Makefile
	cp ${FFLAS_FFPACK_DIR}/include/config-blas.h  /tmp/fflas-ffpack-${VERSION}/include
	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{fflas.h,fflas_fgemm.inl,fflas_fgemv.inl,fflas_ftrsm.inl,fflas_ftrmm.inl,fflas_ftrsv.inl,fflas_fdot.inl,fflas_faxpy.inl,fflas_fcopy.inl,fflas_bounds.inl,fflas_fger.inl,fflas_ftrsm_src.inl,fflas_ftrmm_src.inl,ffpack.h,ffpack_charpoly.inl,ffpack_charpoly_kglu.inl,ffpack_charpoly_danilevski.inl,ffpack_charpoly_kgfast.inl,ffpack_charpoly_kgfastgeneralized.inl,ffpack_frobenius.inl,ffpack_krylovelim.inl,ffpack_ludivine.inl,ffpack_minpoly.inl,unparametric.h,modular-positive.h,modular-balanced.h,modular-int.h}  /tmp/fflas-ffpack-${VERSION}/include/fflas-ffpack

	cp ${FFLAS_FFPACK_DIR}/benchmark/run.sh /tmp/fflas-ffpack-${VERSION}/benchmark
	cp ${FFLAS_FFPACK_DIR}/benchmark/Makefile.Rule.template /tmp/fflas-ffpack-${VERSION}/benchmark/Makefile.Rule
	cp ${FFLAS_FFPACK_DIR}/benchmark/graph/{graph_report.sh,make_graph_file.pl,make_graph.sh} /tmp/fflas-ffpack-${VERSION}/benchmark/graph
	cp ${FFLAS_FFPACK_DIR}/benchmark/html/{fflas.css,html_report.sh,html_report.xsl,process.sh} /tmp/fflas-ffpack-${VERSION}/benchmark/html
	cp ${FFLAS_FFPACK_DIR}/benchmark/src/BLAS_LAPACK/{check-dgemm.C,check-dgetrf.C,check-dgetri.C,check-dtrsm.C,check-dtrtri.C,Makefile} /tmp/fflas-ffpack-${VERSION}/benchmark/src/BLAS_LAPACK
	cp ${FFLAS_FFPACK_DIR}/benchmark/src/FFLAS_FFPACK/{check-fgemm.C,check-ftrsm.C,check-ftrtri.C,check-inverse.C,check-lqup.C,Makefile} /tmp/fflas-ffpack-${VERSION}/benchmark/src/FFLAS_FFPACK
	cp ${FFLAS_FFPACK_DIR}/benchmark/src/BLOCKING/{gnucommand,Makefile,mesure.sh,mulMM.C,plot1-mulMM,tblockmat.C} /tmp/fflas-ffpack-${VERSION}/benchmark/src/BLOCKING
	cp ${FFLAS_FFPACK_DIR}/benchmark/test-src/{mesure-BLAS_LAPACK.sh,mesure-FFLAS_FFPACK.sh,mesure.sh,parameter.in} /tmp/fflas-ffpack-${VERSION}/benchmark/test-src

	tar zcvf fflas-ffpack-${VERSION}.tar.gz  -C /tmp fflas-ffpack-${VERSION}
	rm -rf /tmp/fflas-ffpack-${VERSION}

linbox:

	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{fflas.h,fflas_fgemm.inl,fflas_fgemv.inl,fflas_ftrsm.inl,fflas_ftrsv.inl,fflas_ftrmm.inl,fflas_fdot.inl,fflas_fcopy.inl,fflas_faxpy.inl,fflas_fger.inl,fflas_bounds.inl,fflas_ftrsm_src.inl,fflas_ftrmm_src.inl} ${LINBOX_DIR}/linbox/fflas/
	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{ffpack.h,ffpack_charpoly.inl,ffpack_frobenius.inl,ffpack_ludivine.inl,ffpack_minpoly.inl,ffpack_krylovelim.inl,ffpack_charpoly_danilevski.inl,ffpack_charpoly_kglu.inl,ffpack_charpoly_kgfast.inl,ffpack_charpoly_kgfastgeneralized.inl} ${LINBOX_DIR}/linbox/ffpack/

