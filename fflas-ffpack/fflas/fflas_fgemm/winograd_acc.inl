/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2014 the LinBox group
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

/** @file fflas/fflas_fgemm/winograd_acc.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_acc_INL
#define __FFLASFFPACK_fgemm_winograd_acc_INL

namespace FFLAS { namespace BLAS3 {

	//! @bug this one fails
	template < class Field >
	inline void WinogradAccOld (const Field& F,
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

		FFLASFFPACK_check(!F.isZero(beta));

		typename Field::Element mbeta  ;
		F.neg(mbeta,beta);

		size_t imaxb, jmaxb, imaxa, jmaxa, ldx2;
		size_t x3rd = std::max(mr,kr);
		const typename Field::Element* d11,*d12,*d21,*d22;
		typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2,*dx3;
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C+nr+mr*ldc;



		size_t ldx3;
		// Three temporary submatrices are required
		typename Field::Element* X1 = new typename Field::Element[mr*nr];
		typename Field::Element* X2 = new typename Field::Element[mr*kr];
		typename Field::Element* X3 = new typename Field::Element[x3rd*nr];

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			imaxa = kr;
			ldx2 = jmaxa = mr;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			imaxa = mr;
			ldx2 = jmaxa = kr;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			imaxb = nr;
			jmaxb = kr;
			ldx3 = x3rd;
		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			imaxb = kr;
			ldx3 = jmaxb = nr;
		}

		// P2 = alpha . A12 * B21 + beta . C11  in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, kmax,w-1,base);

		// T3 = B22 - B12 in X3
		d12 = B12; d22 = B22; dx3 = X3;
		fsub(F,imaxb,jmaxb,d22,ldb,d12,ldb,dx3,ldx3);

		// S3 = A11 - A21 in X2
		d11 = A11; d21 = A21; dx2 = X2;
		fsub(F,imaxa,jmaxa,d11,lda,d21,lda,dx2,ldx2);

		// C22 = C22 - C12 if beta != 0
		d12c = C12;
		d22c = C22;
		fsubin(F,mr,nr,d12c,ldc,d22c,ldc);

		// C21 = C21 - C22
		d21c = C21;
		d22c = C22;
		fsubin(F,mr,nr,d22c,ldc,d21c,ldc);

		// P7 = alpha . S3 * T3 + beta . C22 in C22
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, beta, C22, ldc, kmax, w-1,base);

		// T1 = B12 - B11 in X3
		d11 = B11; d12 = B12; dx3 = X3;
		fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx3,ldx3);

		// S1 = A21 + A22 in X2
		d21 = A21; d22 = A22; dx2 = X2;
		fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx2,ldx2);

		// P5 = alpha . S1*T1 + beta . C12 in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, beta, C12, ldc, kmax, w-1,base);

		// T2 = B22 - T1 in X3
		d22 = B22; dx3 = X3;
		fsub(F,imaxb,jmaxb,d22,ldb,dx3,ldx3,dx3,ldx3);

		// S2 = S1 - A11 in X2
		d11 = A11; dx2 = X2;
		fsubin(F,imaxa,jmaxa,d11,lda,dx2,ldx2);

		// P6 = alpha . S2 * T2 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.zero, X1, nr, kmax, w-1,base);

		// T4 = T2 - B21 in X3
		d21 = B21;dx3=X3;
		fsubin(F,imaxb,jmaxb,d21,ldb,dx3,ldx3);

		// S4 = A12 -S2 in X2
		d12 = A12; dx2 = X2;
		fsub(F,imaxa,jmaxa,d12,lda,dx2,ldx2,dx2,ldx2);

		// P4 = alpha . A22 * T4 - beta . C21 in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldx3, mbeta, C21, ldc, kmax, w-1,base);

		// P1 = alpha . A11 * B11 in X3
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X3, nr, kmax, w-1,base);

		//  U1 = P2 + P1 in C11
		d11c = C11; dx3 = X3;
		faddin(F,mr,nr,dx3,nr,d11c,ldc);

		d12c = C12; dx1=X1; dx3=X3; d21c = C21; d22c = C22;

#if 1
		// U2 = P1 + P6 in tmpU2/dx1  and
		faddin(F, mr, nr, dx3, nr, dx1, nr);

		// U3 = P7 + U2 in tmpU3/dx3  and
		fadd(F, mr, nr, dx1, nr, d22c, ldc, dx3, nr);

		// U7 = P5 + U3 in C22    and
		fadd(F, mr, nr, d12c, ldc, dx3, nr, d22c, ldc);

		// U4 = P5 + U2 in C12    and
		faddin(F, mr, nr, dx1, nr, d12c, ldc);

		// U6 = U3 - P4 in C21    and
		fsub(F, mr, nr, dx3, nr, d21c, ldc, d21c, ldc);

#else /*  not using functions */
		for (size_t i = 0; i < mr;
		     ++i, d12c += ldc, dx1 += nr, dx3 += nr, d22c+=ldc, d21c += ldc) {
			for (size_t j=0;j < nr;++j) {
				F.add (tmpU2, *(dx3 + j), *(dx1 + j));    // temporary U2 = P1 + P6
				F.add (tmpU3, *(dx1+j), *(d22c + j));      // temporary U3 = U2 + P7
				F.add (*(d22c + j), *(d12c + j), *(dx3+j));  // U7 = P5 + U3 in C22
				F.addin (*(d12c + j), *(dx1+j));             // U4 = P5 + U2 in C12
				F.sub (*(d21c + j), *(dx3+j), *(d21c + j)); // U6 = U3 - P4 in C21
			}
		}
#endif

		delete[] X1;
		delete[] X3;

		// P3 = alpha . S4*B22 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, B22, ldb, F.one, C12, ldc, kmax, w-1,base);

		// U5 = P3 + U4 in C12
#if 0
		d12c = C12; dx1 = X1;
		for (size_t i = 0; i < mr; ++i, d12c += ldc, dx1 += nr)
			for (size_t j = 0; j < nr; ++j)
				F.addin (*(d12c + j), *(dx1 + j));
#endif
		delete[] X2;

	} // WinogradAccOld

	template < class Field >
	inline void WinogradAcc (const Field& F,
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

		FFLASFFPACK_check(!F.isZero(beta));

		typename Field::Element mbeta  ;
		F.neg(mbeta,beta);

		size_t imaxb, jmaxb, imaxa, jmaxa, ldx2;
		size_t x3rd = std::max(mr,kr);
		const typename Field::Element* d11,*d12,*d21,*d22;
		typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2,*dx3;
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C+nr+mr*ldc;



		size_t ldx3;
		// Three temporary submatrices are required
		typename Field::Element* X1 = new typename Field::Element[mr*nr];
		typename Field::Element* X2 = new typename Field::Element[mr*kr];
		typename Field::Element* X3 = new typename Field::Element[x3rd*nr];

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			imaxa = kr;
			ldx2 = jmaxa = mr;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			imaxa = mr;
			ldx2 = jmaxa = kr;
		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			imaxb = nr;
			jmaxb = kr;
			ldx3 = x3rd;
		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			imaxb = kr;
			ldx3 = jmaxb = nr;
		}


		// T1 = B12 - B11 in X3
		d11 = B11; d12 = B12; dx3 = X3;
		fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx3,ldx3);

		// S1 = A21 + A22 in X2
		d21 = A21; d22 = A22; dx2 = X2;
		fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx2,ldx2);

		// P5 = alpha . S1*T1 + beta . C12 in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.zero, X1, nr, kmax, w-1,base);

		// C22 = P5 + beta C22 in C22
		d22c = C22; dx1 = X1;
		fadd(F,mr,nr,dx1,nr,beta,d22c,ldc,d22c,ldc);

		// C12 = P5 + beta C12 in C12
		dx1 = X1; d12c = C12;
		fadd(F,mr,nr,dx1,nr,beta,d12c,ldc,d12c,ldc);

		// P1 = alpha . A11 * B11 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, kmax, w-1,base);


		// P2 = alpha . A12 * B21 + beta . C11  in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, kmax,w-1,base);

		//  U1 = P2 + P1 in C11
		d11c = C11; dx1 = X1;
		faddin(F,mr,nr,dx1,nr,d11c,ldc);

		// T2 = B22 - T1 in X3
		d22 = B22; dx3 = X3;
		fsub(F,imaxb,jmaxb,d22,ldb,dx3,ldx3,dx3,ldx3);

		// S2 = S1 - A11 in X2
		d11 = A11; dx2 = X2;
		fsubin(F,imaxa,jmaxa,d11,lda,dx2,ldx2);

		// U2 = P6 + P1 = alpha . S2 * T2 + P1 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.one, X1, nr, kmax, w-1,base);

		// U4 = U2 + P5 in C12
		d12c = C12; dx1 = X1;
		faddin(F,mr,nr,dx1,nr,d12c,ldc);

		// T4 = T2 - B21 in X3
		d21 = B21;dx3=X3;
		fsubin(F,imaxb,jmaxb,d21,ldb,dx3,ldx3);

		// S4 = A12 -S2 in X2
		d12 = A12; dx2 = X2;
		fsub(F,imaxa,jmaxa,d12,lda,dx2,ldx2,dx2,ldx2);

		// P4 = alpha . A22 * T4 - beta . C21 in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldx3, mbeta, C21, ldc, kmax, w-1,base);

		// U5 = P3 + U4 = alpha . S4*B22 + U4 in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, B22, ldb, F.one, C12, ldc, kmax, w-1,base);

		// T3 = B22 - B12 in X3
		d12 = B12; d22 = B22; dx3 = X3;
		fsub(F,imaxb,jmaxb,d22,ldb,d12,ldb,dx3,ldx3);

		// S3 = A11 - A21 in X2
		d11 = A11; d21 = A21; dx2 = X2;
		fsub(F,imaxa,jmaxa,d11,lda,d21,lda,dx2,ldx2);

		// U3 = P7 + U2  = alpha . S3 * T3 + U2 in X1
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.one, X1, nr, kmax, w-1,base);

		// U7 =  U3 + C22 in C22
		d22c = C22; dx1 = X1; d12c = C12;
		faddin(F,mr,nr,dx1,nr,d22c,ldc);

		// U6 = U3 - P4 in C21
		dx1 = X1; d21c = C21;
		fsub(F,mr,nr,dx1,nr,d21c,ldc,d21c,ldc);
		delete[] X2;

	} // WinogradAcc


} // BLAS3

} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_acc_INL

