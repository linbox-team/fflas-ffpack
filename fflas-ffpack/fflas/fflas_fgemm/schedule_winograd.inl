/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            BB <bbboyer@ncsu.edu>
 *
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

	template < class Field >
	inline void Winograd (const Field& F,
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element alpha,
			      const typename Field::Element* A,const size_t lda,
			      const typename Field::Element* B,const size_t ldb,
			      const typename Field::Element  beta,
			      typename Field::Element * C, const size_t ldc,
			      // const size_t kmax, const size_t w, const FFLAS_BASE base
			      const Winograd2Helper & WH
			      )
	{

		FFLASFFPACK_check(F.isZero(beta));

		size_t lb, cb, la, ca, ldX2;
		// size_t x3rd = std::max(mr,kr);
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C21+nr;


		size_t x1rd = std::max(nr,kr);
		size_t ldX1;
		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			la = kr;
			ca = mr;
			ldX1 = mr;
		}
		else {
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
		}
		else {
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			lb = kr;
			ldX2 = cb = nr;
		}


		Winograd2Helper H = WH ;
		H.w = H.w -1 ;
		// Two temporary submatrices are required
		typename Field::Element* X2 = new typename Field::Element[kr*nr];

		// T3 = B22 - B12 in X2
		fsub(F,lb,cb,B22,ldb,B12,ldb,X2,ldX2);

		// S3 = A11 - A21 in X1
		typename Field::Element* X1 = new typename Field::Element[mr*x1rd];
		fsub(F,la,ca,A11,lda,A21,lda,X1,ldX1);

		// P7 = alpha . S3 * T3  in C21
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C21, ldc, H);

		// T1 = B12 - B11 in X2
		fsub(F,lb,cb,B12,ldb,B11,ldb,X2,ldX2);

		// S1 = A21 + A22 in X1

		fadd(F,la,ca,A21,lda,A22,lda,X1,ldX1);

		// P5 = alpha . S1*T1 in C22
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C22, ldc, H);

		// T2 = B22 - T1 in X2
		fsub(F,lb,cb,B22,ldb,X2,ldX2,X2,ldX2);

		// S2 = S1 - A11 in X1
		fsubin(F,la,ca,A11,lda,X1,ldX1);

		// P6 = alpha . S2 * T2 in C12
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, X2, ldX2, F.zero, C12, ldc, H);

		// S4 = A12 -S2 in X1
		fsub(F,la,ca,A12,lda,X1,ldX1,X1,ldX1);

		// P3 = alpha . S4*B22 in C11
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, X1, ldX1, B22, ldb, F.zero, C11, ldc, H);

		// P1 = alpha . A11 * B11 in X1
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, H);

		// U2 = P1 + P6 in tmpU2  and
		faddin(F,mr,nr,X1,nr,C12,ldc);

		// U3 = P7 + U2 in tmpU3  and
		faddin(F,mr,nr,C12,ldc,C21,ldc);

		// U7 = P5 + U3 in C22    and
		faddin(F,mr,nr,C22,ldc,C12,ldc);

		// U4 = P5 + U2 in C12    and
		faddin(F,mr,nr,C21,ldc,C22,ldc);

		// U5 = P3 + U4 in C12
		faddin(F,mr,nr,C11,ldc,C12,ldc);

		// T4 = T2 - B21 in X2
		fsubin(F,lb,cb,B21,ldb,X2,ldX2);

		// P4 = alpha . A22 * T4 in C11
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, A22, lda, X2, ldX2, F.zero, C11, ldc, H);

		delete[] X2;
		// U6 = U3 - P4 in C21
		fsubin(F,mr,nr,C11,ldc,C21,ldc);

		// P2 = alpha . A12 * B21  in C11
		fgemm2 (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, H);

		//  U1 = P2 + P1 in C11
		faddin(F,mr,nr,X1,nr,C11,ldc);

		delete[] X1;

	} // Winograd

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_INL

