/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Written by BB <bbboyer@ncsu.edu>
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
			      const size_t kmax, const size_t w, const FFLAS_BASE base)
	{

		FFLASFFPACK_check(F.isZero(beta));

		size_t imaxb, jmaxb, imaxa, jmaxa, ldx2;
		// size_t x3rd = std::max(mr,kr);
		const typename Field::Element* d11,*d12,*d21,*d22;
		typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2;
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C+nr+mr*ldc;


		size_t x1rd = std::max(nr,kr);
		size_t ldx1;
		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			imaxa = kr;
			jmaxa = mr;
			ldx1 = mr;
		}
		else {
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			imaxa = mr;
			jmaxa = kr;
			ldx1  = x1rd;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			imaxb = nr;
			jmaxb = kr;
			ldx2 = kr;
		}
		else {
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			imaxb = kr;
			ldx2 = jmaxb = nr;
		}


		// Two temporary submatrices are required

		typename Field::Element* X2 = new typename Field::Element[kr*nr];

		// T3 = B22 - B12 in X2
		d12 = B12; d22 = B22; dx2 = X2;
		fsub(F,imaxb,jmaxb,d22,ldb,d12,ldb,dx2,ldx2);

		// S3 = A11 - A21 in X1
		typename Field::Element* X1 = new typename Field::Element[mr*x1rd];
		d11 = A11; d21 = A21; dx1 = X1;
		fsub(F,imaxa,jmaxa,d11,lda,d21,lda,dx1,ldx1);

		// P7 = alpha . S3 * T3  in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C21, ldc, kmax, w-1, base);

		// T1 = B12 - B11 in X2
		d11 = B11; d12 = B12; dx2 = X2;
		fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx2,ldx2);

		// S1 = A21 + A22 in X1

		d21 = A21; d22 = A22; dx1 = X1;
		fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx1,ldx1);

		// P5 = alpha . S1*T1 in C22
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C22, ldc, kmax, w-1, base);

		// T2 = B22 - T1 in X2
		d22 = B22; dx2 = X2;
		fsub(F,imaxb,jmaxb,d22,ldb,dx2,ldx2,dx2,ldx2);

		// S2 = S1 - A11 in X1
		d11 = A11; dx1 = X1;
		fsubin(F,imaxa,jmaxa,d11,lda,dx1,ldx1);

		// P6 = alpha . S2 * T2 in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C12, ldc, kmax, w-1, base);

		// S4 = A12 -S2 in X1
		d12 = A12; dx1 = X1;
		fsub(F,imaxa,jmaxa,d12,lda,dx1,ldx1,dx1,ldx1);

		// P3 = alpha . S4*B22 in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, B22, ldb, F.zero, C11, ldc, kmax, w-1, base);

		// P1 = alpha . A11 * B11 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, kmax, w-1, base);



		// U2 = P1 + P6 in tmpU2  and
		// U3 = P7 + U2 in tmpU3  and
		// U7 = P5 + U3 in C22    and
		// U4 = P5 + U2 in C12    and
		d12c = C12; dx1=X1; d21c = C21; d22c = C22;
		for (size_t i = 0; i < mr;
		     ++i, d12c += ldc, dx1 += nr, d22c+=ldc, d21c += ldc) {
			for (size_t j=0;j < nr;++j) {
				F.addin ( *(d12c + j), *(dx1 +j));    // U2 = P1 + P6
				F.addin ( *(d21c+j)  , *(d12c+j));      //  U3 = U2 + P7
				F.addin ( *(d12c + j), *(d22c+j));   // U4 = P5 + U2 in C12
				F.addin ( *(d22c + j), *(d21c+j));  // U7 = P5 + U3 in C22
			}
		}

		// U5 = P3 + U4 in C12
		d12c = C12; d11 = C11;
		faddin(F,mr,nr,d11,ldc,d12c,ldc);

		// T4 = T2 - B21 in X2
		d21 = B21;dx2=X2;
		fsubin(F,imaxb,jmaxb,d21,ldb,dx2,ldx2);

		// P4 = alpha . A22 * T4 in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X2, ldx2, F.zero, C11, ldc, kmax, w-1, base);

		delete[] X2;
		// U6 = U3 - P4 in C21
		d21c = C21; d11c = C11;
		fsubin(F,mr,nr,d11c,ldc,d21c,ldc);

		// P2 = alpha . A12 * B21  in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, kmax,w-1, base);

		//  U1 = P2 + P1 in C11
		d11c = C11; dx1 = X1;
		faddin(F,mr,nr,dx1,nr,d11c,ldc);

		delete[] X1;

	} // Winograd

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_INL

