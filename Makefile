VERSION=1.3.0

FFLAS_FFPACK_DIR=/home/pernet/Logiciels/fflas-ffpack
LINBOX_DIR=/home/pernet/Logiciels/linbox


all:
	cd tests
	make

dist:
	mkdir /tmp/fflas-ffpack-${VERSION}
	mkdir /tmp/fflas-ffpack-${VERSION}/tests /tmp/fflas-ffpack-${VERSION}/include /tmp/fflas-ffpack-${VERSION}/include/fflas-ffpack
	cp ${FFLAS_FFPACK_DIR}/{AUTHORS,ChangeLog,COPYING,INSTALL,README} /tmp/fflas-ffpack-${VERSION}/
	cp ${FFLAS_FFPACK_DIR}/tests/{timer.h,timer.C,Matio.h,test-fgemm.C,test-fgemv.C,test-lqup.C,test-charpoly.C,test-compressQ.C,test-frobenius.C,test-fsquare.C,test-det.C,test-invert.C,test-krylov-elim.C,test-rank.C,testeur_fgemm.C,testeur_lqup.C,testeur_ftrsm.C,dense_generator.C}  /tmp/fflas-ffpack-${VERSION}/tests/
	cp ${FFLAS_FFPACK_DIR}/tests/Makefile.template /tmp/fflas-ffpack-${VERSION}/tests/Makefile
	cp ${FFLAS_FFPACK_DIR}/include/config-blas.h  /tmp/fflas-ffpack-${VERSION}/include
	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{fflas.h,fflas_fgemm.inl,fflas_fgemv.inl,fflas_ftrsm.inl,fflas_ftrmm.inl,fflas_ftrsv.inl,fflas_fdot.inl,fflas_faxpy.inl,fflas_fcopy.inl,fflas_bounds.inl,fflas_fger.inl,fflas_ftrsm_src.inl,fflas_ftrmm_src.inl,ffpack.h,ffpack_charpoly.inl,ffpack_charpoly_kglu.inl,ffpack_charpoly_danilevski.inl,ffpack_charpoly_kgfast.inl,ffpack_charpoly_kgfastgeneralized.inl,ffpack_frobenius.inl,ffpack_krylovelim.inl,ffpack_ludivine.inl,ffpack_minpoly.inl,unparametric.h,modular-positive.h,modular-balanced.h,modular-int.h}  /tmp/fflas-ffpack-${VERSION}/include/fflas-ffpack
	tar zcvf fflas-ffpack-${VERSION}.tar.gz  -C /tmp fflas-ffpack-${VERSION}
	rm -rf /tmp/fflas-ffpack-${VERSION}

linbox:

	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{fflas.h,fflas_fgemm.inl,fflas_fgemv.inl,fflas_ftrsm.inl,fflas_ftrsv.inl,fflas_ftrmm.inl,fflas_fdot.inl,fflas_fcopy.inl,fflas_faxpy.inl,fflas_fger.inl,fflas_bounds.inl,fflas_ftrsm_src.inl,fflas_ftrmm_src.inl} ${LINBOX_DIR}/linbox/fflas/
	cp ${FFLAS_FFPACK_DIR}/include/fflas-ffpack/{ffpack.h,ffpack_charpoly.inl,ffpack_frobenius.inl,ffpack_ludivine.inl,ffpack_minpoly.inl,ffpack_krylovelim.inl,ffpack_charpoly_danilevski.inl,ffpack_charpoly_kglu.inl,ffpack_charpoly_kgfast.inl,ffpack_charpoly_kgfastgeneralized.inl} ${LINBOX_DIR}/linbox/ffpack/

