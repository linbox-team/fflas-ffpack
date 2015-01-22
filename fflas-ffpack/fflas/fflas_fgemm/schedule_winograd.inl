/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            BB <bbboyer@ncsu.edu>
 *            Ziad Sultan <ziad.sultan@imag.fr>
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas/fflas_fgemm/winograd.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_INL
#define __FFLASFFPACK_fgemm_winograd_INL

namespace FFLAS { namespace BLAS3 {

		template < class Field, class FieldTrait >
	inline void WinoPar (const Field& F,
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element alpha,
			      typename Field::ConstElement_ptr A,const size_t lda,
			      typename Field::ConstElement_ptr B,const size_t ldb,
			      const typename Field::Element  beta,
			      typename Field::Element_ptr C, const size_t ldc,
			      // const size_t kmax, const size_t w, const FFLAS_BASE base
			      MMHelper<Field, MMHelperAlgo::WinogradPar, FieldTrait> & WH
			     )
	{

		FFLASFFPACK_check(F.isZero(beta));

		typedef MMHelper<Field, MMHelperAlgo::WinogradPar, FieldTrait > MMH_t;

		const typename MMH_t::DelayedField & DF = WH.delayedField;

		size_t lb, cb, la, ca, ldX2;
		    // size_t x3rd = std::max(mr,kr);
		typename Field::ConstElement_ptr A11=A, A12, A21, A22;
		typename Field::ConstElement_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

		size_t x1rd = std::max(nr,kr);
		size_t ldX1;
		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr;
			ca = mr;
			ldX1 = mr;
		} else {
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr;
			ca = kr;
			ldX1  = x1rd;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr;
			cb = kr;
			ldX2 = kr;
		} else {
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr;
			ldX2 = cb = nr;
		}

		//		std::cout<<"mr: "<<mr<<" nr: "<<nr<<" kr: "<<kr<<" x1rd: "<<x1rd<<std::endl;
		// 11 temporary submatrices are required
		typename Field::Element_ptr X21 = fflas_new (F, kr, nr);
		typename Field::Element_ptr X11 = fflas_new (F,mr,x1rd);
		
		typename Field::Element_ptr X22 = fflas_new (F, kr, nr);
		typename Field::Element_ptr X12 = fflas_new (F,mr,x1rd);

		typename Field::Element_ptr X23 = fflas_new (F, kr, nr);
		typename Field::Element_ptr X13 = fflas_new (F,mr,x1rd);

		typename Field::Element_ptr X24 = fflas_new (F, kr, nr);
		typename Field::Element_ptr X14 = fflas_new (F,mr,x1rd);
		typename Field::Element_ptr X15 = fflas_new (F,mr,x1rd);

		typename Field::Element_ptr C_11 = fflas_new (F,mr,nr);
		typename Field::Element_ptr CC_11 = fflas_new (F,mr,nr);


		// P1 = alpha . A11 * B11 in X1
		MMH_t H1(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
		TASK(MODE(READ(A11, B11) WRITE(X15) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X15, x1rd, H1););
		
		// T3 = B22 - B12 in X21  and S3 = A11 - A21 in X11
		TASK(MODE(READ(B22, B12) WRITE(X21) REFERENCE(DF)),
		     fsub(DF,lb,cb,B22,ldb,B12,ldb,X21,ldX2););
		TASK(MODE(READ(A11, A21) WRITE(X11) REFERENCE(DF)),
		     fsub(DF,la,ca,A11,lda,A21,lda,X11,ldX1););

		// T1 = B12 - B11 in X22 and  S1 = A21 + A22 in X12
		TASK(MODE(READ(B11, B12) WRITE(X22) REFERENCE(DF)),
		     fsub(DF,lb,cb,B12,ldb,B11,ldb,X22,ldX2););
		TASK(MODE(READ(A12, A22) WRITE(X12) REFERENCE(DF)),
		     fadd(DF,la,ca,A21,lda,A22,lda,X12,ldX1););

		CHECK_DEPENDENCIES;

		// T2 = B22 - T1 in X23 and  S2 = S1 - A11 in X13
		TASK(MODE(READ(B22, X22) READWRITE(X23) REFERENCE(DF)),
		     fsub(DF,lb,cb,B22,ldb,X22,ldX2,X23,ldX2););
		TASK(MODE(READ(A11, X12) READWRITE(X13) REFERENCE(DF)),
		     //		     fsub(DF,la,ca,A11,lda,X12,ldX1,X13,ldX1););
		     fsub(DF,la,ca,X12,ldX1,A11,lda,X13,ldX1););
		/*
		  fsub(DF,lb,cb,B22,ldb,X2,ldX2,X2,ldX2);
		  fsubin(DF,la,ca,A11,lda,X1,ldX1););
		*/
		CHECK_DEPENDENCIES;

		// T4 = T2 - B21 in X2 and S4 = A12 -S2 in X1
		TASK(MODE(READ(B21, X23) READWRITE(X24) REFERENCE(DF)),
		     //		     fsub(DF,lb,cb,B21,ldb,X23,ldX2,X24,ldX2);
		     fsub(DF,lb,cb,X23,ldX2,B21,ldb,X24,ldX2););
		TASK(MODE(READ(A12, X13) READWRITE(X14) REFERENCE(DF)),
		     fsub(DF,la,ca,A12,lda,X13,ldX1,X14,ldX1););

		     /*
		     fsubin(DF,lb,cb,B21,ldb,X2,ldX2);
		     fsub(DF,la,ca,A12,lda,X1,ldX1,X1,ldX1););
		     */
		CHECK_DEPENDENCIES;
		
		// P7 = alpha . S3 * T3  in C21
		MMH_t H7(F, WH.recLevel-1, -(WH.Amax-WH.Amin), WH.Amax-WH.Amin, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0,0);
		TASK(MODE(READ(X11, X21) WRITE(C21) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, X11, ldX1, X21, ldX2, F.zero, C21, ldc, H7););

		// P5 = alpha . S1*T1 in C22
		MMH_t H5(F, WH.recLevel-1, 2*WH.Amin, 2*WH.Amax, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0, 0);
		TASK(MODE(READ(X12, X22) WRITE(C22) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, X12, ldX1, X22, ldX2, F.zero, C22, ldc, H5););

		// P6 = alpha . S2 * T2 in C12
		MMH_t H6(F, WH.recLevel-1, 2*WH.Amin-WH.Amax, 2*WH.Amax-WH.Amin, 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin, 0, 0);
		TASK(MODE(READ(X13, X23) WRITE(C12) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, X13, ldX1, X23, ldX2, F.zero, C12, ldc, H6););

		// P3 = alpha . S4*B22 in CC_11
		MMH_t H3(F, WH.recLevel-1, 2*WH.Amin-2*WH.Amax, 2*WH.Amax-2*WH.Amin, WH.Bmin, WH.Bmax, 0, 0);
		TASK(MODE(READ(X14, B22) WRITE(CC_11) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, X14, ldX1, B22, ldb, F.zero, CC_11, nr, H3););

		// P4 = alpha . A22 * T4 in C_11
		MMH_t H4(F, WH.recLevel-1, WH.Amin, WH.Amax, 2*WH.Bmin-2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin, 0, 0);
		TASK(MODE(READ(A22) WRITE(C_11) READWRITE(X24, X22, X23, X21) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X24, ldX2, F.zero, C_11, nr, H4);
		     );

				// P2 = alpha . A12 * B21  in C11
		MMH_t H2(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
		TASK(MODE(READ(A12, B21) WRITE(C11) REFERENCE(F)),
		     fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, H2););
		CHECK_DEPENDENCIES;

		double U2Min, U2Max;
		double U3Min, U3Max;
		double U4Min, U4Max;
		double U7Min, U7Max;
		double U5Min, U5Max;
		// U2 = P1 + P6 in C12  and
		// U3 = P7 + U2 in C21  and
		// U4 = P5 + U2 in C12    and
		// U7 = P5 + U3 in C22    and
		// U5 = P3 + U4 in C12
		// BIG TASK with 5 Addin function calls
		TASK(MODE(READWRITE(X15, C12) REFERENCE(DF, F, WH)),
		     if (Protected::NeedPreAddReduction(U2Min, U2Max, H1.Outmin, H1.Outmax, H6.Outmin, H6.Outmax, WH)){
			     freduce (F, mr, x1rd, X15, x1rd);
			     freduce (F, mr, nr, C12, ldc);
		     }
		     faddin(DF,mr,nr,X15,x1rd,C12,ldc);
		     );
		CHECK_DEPENDENCIES;
		TASK(MODE(READWRITE(C12, C21) REFERENCE(DF, F, WH)),
		     if (Protected::NeedPreAddReduction(U3Min, U3Max, U2Min, U2Max, H7.Outmin, H7.Outmax, WH)){
			     freduce (F, mr, nr, C12, ldc);
			     freduce (F, mr, nr, C21, ldc);
		     }
		     faddin(DF,mr,nr,C12,ldc,C21,ldc);
		     );
		CHECK_DEPENDENCIES;
		TASK(MODE(READWRITE(C12, C22) REFERENCE(DF, F, WH)),
		     if (Protected::NeedPreAddReduction(U4Min, U4Max, U2Min, U2Max, H5.Outmin, H5.Outmax, WH)){
			     freduce (F, mr, nr, C22, ldc);
			     freduce (F, mr, nr, C12, ldc);
		     }
		     faddin(DF,mr,nr,C22,ldc,C12,ldc);
		     );
		CHECK_DEPENDENCIES;
		TASK(MODE(READWRITE(C22, C21) REFERENCE(DF, F, WH)),
		     if (Protected::NeedPreAddReduction (U7Min,U7Max, U3Min, U3Max, H5.Outmin,H5.Outmax, WH) ){
			     freduce (F, mr, nr, C21, ldc);
			     freduce (F, mr, nr, C22, ldc);
		     }
		     faddin(DF,mr,nr,C21,ldc,C22,ldc);
		     );
		TASK(MODE(READWRITE(C12, CC_11) REFERENCE(DF, F, WH)),
		     if (Protected::NeedPreAddReduction (U5Min,U5Max, U4Min, U4Max, H3.Outmin, H3.Outmax, WH) ){
			     freduce (F, mr, nr, C12, ldc);
			     freduce (F, mr, nr, CC_11, nr);
		     }
		     faddin(DF,mr,nr,CC_11,nr,C12,ldc);
		     );
		CHECK_DEPENDENCIES;

		// U6 = U3 - P4 in C21
		double U6Min, U6Max;
		TASK(MODE(READWRITE(C_11, C21) REFERENCE(DF, F, WH)),
		if (Protected::NeedPreSubReduction (U6Min,U6Max, U3Min, U3Max, H4.Outmin,H4.Outmax, WH) ){
			freduce (F, mr, nr, C_11, nr);
			freduce (F, mr, nr, C21, ldc);
		}
		     fsubin(DF,mr,nr,C_11,nr,C21,ldc););
		
		    //CHECK_DEPENDENCIES;

		//  U1 = P2 + P1 in C11
		double U1Min, U1Max;
		TASK(MODE(READWRITE(C11, X15/*, X14, X13, X12, X11*/) REFERENCE(DF, F, WH)),
		if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.Outmin, H1.Outmax, H2.Outmin,H2.Outmax, WH) ){
			freduce (F, mr, nr, X15, x1rd);
			freduce (F, mr, nr, C11, ldc);
		}
		     faddin(DF,mr,nr,X15,x1rd,C11,ldc);
		
		     );

		WH.Outmin = std::min (U1Min, std::min (U5Min, std::min (U6Min, U7Min)));
		WH.Outmax = std::max (U1Max, std::max (U5Max, std::max (U6Max, U7Max)));

		WAIT;
		
		fflas_delete (CC_11);
 		fflas_delete (C_11);
		fflas_delete (X15);
		fflas_delete (X14);
		fflas_delete (X24);
		fflas_delete (X13);
		fflas_delete (X23);
		fflas_delete (X12);
		fflas_delete (X22);
		fflas_delete (X11);
		fflas_delete (X21);
		     

	} //wino parallel


		template < class Field, class FieldTrait >
		inline void Winograd (const Field& F,
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element alpha,
			      typename Field::ConstElement_ptr A,const size_t lda,
			      typename Field::ConstElement_ptr B,const size_t ldb,
			      const typename Field::Element  beta,
			      typename Field::Element_ptr C, const size_t ldc,
			      // const size_t kmax, const size_t w, const FFLAS_BASE base
			      MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait> & WH
			     )
	{

		FFLASFFPACK_check(F.isZero(beta));

		typedef MMHelper<Field, MMHelperAlgo::Winograd, FieldTrait > MMH_t;

		const typename MMH_t::DelayedField & DF = WH.delayedField;

		size_t lb, cb, la, ca, ldX2;
		    // size_t x3rd = std::max(mr,kr);
		typename Field::ConstElement_ptr A11=A, A12, A21, A22;
		typename Field::ConstElement_ptr B11=B, B12, B21, B22;
		typename Field::Element_ptr C11=C, C12=C+nr, C21=C+mr*ldc, C22=C21+nr;

		size_t x1rd = std::max(nr,kr);
		size_t ldX1;
		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr;
			ca = mr;
			ldX1 = mr;
		} else {
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			la = mr;
			ca = kr;
			ldX1  = x1rd;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			lb = nr;
			cb = kr;
			ldX2 = kr;
		} else {
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr;
			ldX2 = cb = nr;
		}
		//		std::cout<<"mr: "<<mr<<" nr: "<<nr<<" kr: "<<kr<<" x1rd: "<<x1rd<<" ldX1: "<<ldX1<<" ldX2: "<<ldX2<<std::endl;
		// Two temporary submatrices are required
		typename Field::Element_ptr X2 = fflas_new (F, kr, nr);

		// T3 = B22 - B12 in X2
		fsub(DF,lb,cb,B22,ldb,B12,ldb,X2,ldX2);

		// S3 = A11 - A21 in X1
		typename Field::Element_ptr X1 = fflas_new (F,mr,x1rd);
		fsub(DF,la,ca,A11,lda,A21,lda,X1,ldX1);

		// P7 = alpha . S3 * T3  in C21
		MMH_t H7(F, WH.recLevel-1, -(WH.Amax-WH.Amin), WH.Amax-WH.Amin, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0,0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C21, ldc, H7);

		// T1 = B12 - B11 in X2
		fsub(DF,lb,cb,B12,ldb,B11,ldb,X2,ldX2);

		// S1 = A21 + A22 in X1

		fadd(DF,la,ca,A21,lda,A22,lda,X1,ldX1);

		// P5 = alpha . S1*T1 in C22
		MMH_t H5(F, WH.recLevel-1, 2*WH.Amin, 2*WH.Amax, -(WH.Bmax-WH.Bmin), WH.Bmax-WH.Bmin, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C22, ldc, H5);

		// T2 = B22 - T1 in X2
		fsub(DF,lb,cb,B22,ldb,X2,ldX2,X2,ldX2);

		// S2 = S1 - A11 in X1
		fsubin(DF,la,ca,A11,lda,X1,ldX1);

		// P6 = alpha . S2 * T2 in C12
		MMH_t H6(F, WH.recLevel-1, 2*WH.Amin-WH.Amax, 2*WH.Amax-WH.Amin, 2*WH.Bmin-WH.Bmax, 2*WH.Bmax-WH.Bmin, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C12, ldc, H6);

		// S4 = A12 -S2 in X1
		fsub(DF,la,ca,A12,lda,X1,ldX1,X1,ldX1);

		// P3 = alpha . S4*B22 in C11
		MMH_t H3(F, WH.recLevel-1, 2*WH.Amin-2*WH.Amax, 2*WH.Amax-2*WH.Amin, WH.Bmin, WH.Bmax, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, B22, ldb, F.zero, C11, ldc, H3);

		// P1 = alpha . A11 * B11 in X1
		MMH_t H1(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, H1);

		// U2 = P1 + P6 in C12  and
		double U2Min, U2Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction(U2Min, U2Max, H1.Outmin, H1.Outmax, H6.Outmin, H6.Outmax, WH)){
			freduce (F, mr, nr, X1, nr);
			freduce (F, mr, nr, C12, ldc);
		}
		faddin(DF,mr,nr,X1,nr,C12,ldc);

		// U3 = P7 + U2 in C21  and
		double U3Min, U3Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction(U3Min, U3Max, U2Min, U2Max, H7.Outmin, H7.Outmax, WH)){
			freduce (F, mr, nr, C12, ldc);
			freduce (F, mr, nr, C21, ldc);
		}
		faddin(DF,mr,nr,C12,ldc,C21,ldc);


		// U4 = P5 + U2 in C12    and
		double U4Min, U4Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction(U4Min, U4Max, U2Min, U2Max, H5.Outmin, H5.Outmax, WH)){
			freduce (F, mr, nr, C22, ldc);
			freduce (F, mr, nr, C12, ldc);
		}
		faddin(DF,mr,nr,C22,ldc,C12,ldc);


		// U7 = P5 + U3 in C22    and
		double U7Min, U7Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction (U7Min,U7Max, U3Min, U3Max, H5.Outmin,H5.Outmax, WH) ){
			freduce (F, mr, nr, C21, ldc);
			freduce (F, mr, nr, C22, ldc);
		}
		faddin(DF,mr,nr,C21,ldc,C22,ldc);


		// U5 = P3 + U4 in C12
		double U5Min, U5Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction (U5Min,U5Max, U4Min, U4Max, H3.Outmin, H3.Outmax, WH) ){
			freduce (F, mr, nr, C12, ldc);
			freduce (F, mr, nr, C11, ldc);
		}
		faddin(DF,mr,nr,C11,ldc,C12,ldc);

		// T4 = T2 - B21 in X2
		fsubin(DF,lb,cb,B21,ldb,X2,ldX2);

		// P4 = alpha . A22 * T4 in C11
		MMH_t H4(F, WH.recLevel-1, WH.Amin, WH.Amax, 2*WH.Bmin-2*WH.Bmax, 2*WH.Bmax-2*WH.Bmin, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A22, lda, X2, ldX2, F.zero, C11, ldc, H4);

		fflas_delete (X2);

		// U6 = U3 - P4 in C21
		double U6Min, U6Max;
		    // This test will be optimized out
		if (Protected::NeedPreSubReduction (U6Min,U6Max, U3Min, U3Max, H4.Outmin,H4.Outmax, WH) ){
			freduce (F, mr, nr, C11, ldc);
			freduce (F, mr, nr, C21, ldc);
		}
		fsubin(DF,mr,nr,C11,ldc,C21,ldc);

		// P2 = alpha . A12 * B21  in C11
		MMH_t H2(F, WH.recLevel-1, WH.Amin, WH.Amax, WH.Bmin, WH.Bmax, 0, 0);
		fgemm (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, H2);

		//  U1 = P2 + P1 in C11
		double U1Min, U1Max;
		    // This test will be optimized out
		if (Protected::NeedPreAddReduction (U1Min, U1Max, H1.Outmin, H1.Outmax, H2.Outmin,H2.Outmax, WH) ){
			freduce (F, mr, nr, X1, nr);
			freduce (F, mr, nr, C11, ldc);
		}
		faddin(DF,mr,nr,X1,nr,C11,ldc);

		fflas_delete (X1);



		WH.Outmin = std::min (U1Min, std::min (U5Min, std::min (U6Min, U7Min)));
		WH.Outmax = std::max (U1Max, std::max (U5Max, std::max (U6Max, U7Max)));
		//#endif

	} // Winograd

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_INL

