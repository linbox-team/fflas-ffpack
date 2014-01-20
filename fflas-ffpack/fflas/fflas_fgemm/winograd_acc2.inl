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

/** @file fflas/fflas_fgemm/winograd_acc2.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_acc2_INL
#define __FFLASFFPACK_fgemm_winograd_acc2_INL

namespace FFLAS { namespace BLAS3 {

	// accw
	template < class Field >
	inline void WinogradAcc2 (const Field& F,
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

		typename Field::Element malpha ;
		F.neg(malpha,alpha);

		// A, B and c submatrices
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C21+nr;



		size_t ldxa;           // ld of X when X=Asub
		size_t ldxc = nr;      // ld of X when X=Csub
		size_t ldy ;           // ld of Y
		size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

		// Three temporary submatrices are required
		typename Field::Element* X = new typename Field::Element[mr*std::max(nr,kr)];
		typename Field::Element* Y = new typename Field::Element[nr*kr];

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			ldxa = mr ;
			la = kr ;
			ca = mr ;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			ldxa = kr ;
			la = mr ;
			ca = kr ;

		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			ldy = kr ;
			lb = nr ;
			cb = kr ;

		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			ldy = nr ;
			lb = kr ;
			cb = nr ;
		}

		// Z1 = C22 - C12         in C22
		fsubin(F,mr,nr,C12,ldc,C22,ldc);
		// Z3 = C12-C21           in C12
		fsubin(F,mr,nr,C21,ldc,C12,ldc);
		// S1 = A21 + A22         in X
		fadd(F,la,ca,A21,lda,A22,lda,X,ldxa);
		// T1 = B12 - B11         in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,ldy);
		// P5 = a S1 T1 + b Z3    in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, beta, C12, ldc, kmax,w-1,base);
		// S2 = S1 - A11          in X
		fsubin(F,la,ca,A11,lda,X,ldxa);
		// T2 = B22 - T1          in Y
		fsub(F,lb,cb,B22,ldb,Y,ldy,Y,ldy);
		// P6 = a S2 T2 + b C21   in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, beta, C21, ldc, kmax,w-1,base);
		// S4 = A12 - S2          in X
		fsub(F,la,ca,A12,lda,X,ldxa,X,ldxa);
		// W1 = P5 + beta Z1      in C22
		fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
		// P3 = a S4 B22 + P5     in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, B22, ldb, F.one, C12, ldc, kmax,w-1,base);
		// P1 = a A11 B11         in X
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, ldxc, kmax,w-1,base);
		// U2 = P6 + P1           in C21
		faddin(F,mr,nr,X,ldxc,C21,ldc);
		// P2 = a A12 B21 + b C11 in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, kmax,w-1,base);
		// U1 = P1 + P2           in C11
		faddin(F,mr,nr,X,ldxc,C11,ldc);
		// U5 = U2 + P3           in C12
		faddin(F,mr,nr,C21,ldc,C12,ldc);
		// S3 =  A11 - A21        in X ;
		fsub(F,la,ca,A11,lda,A21,lda,X,ldxa);
		// T3 = B22 - B12         in Y
		fsub(F,lb,cb,B22,ldb,B12,ldb,Y,ldy);
		// U3 = a S3 T3 + U2      in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, F.one, C21, ldc, kmax,w-1,base);
		delete[] X;
		// U7 = U3 + W1           in C22
		faddin(F,mr,nr,C21,ldc,C22,ldc);
		// T1_ = B12 - B11        in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,ldy);
		// T2_ = B22 - T1_        in Y
		fsub(F,lb,cb,B22,ldb,Y,ldy,Y,ldy);
		// T4 = T2_ - B21         in Y
		fsub(F,lb,cb,Y,ldy,B21,ldb,Y,ldy);
		// U6 = -a A22 T4 + U3    in C21;
		Protected::WinoMain (F, ta, tb, mr, nr, kr, malpha, A22, lda, Y, ldy, F.one, C21, ldc, kmax,w-1,base);
		delete[] Y;


	} // WinogradAccOld

	// accw2
	template < class Field >
	inline void WinogradAcc3 (const Field& F,
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

		typename Field::Element malpha ;
		F.neg(malpha,alpha);

		// A, B and c submatrices
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C21+nr;



		size_t ldxa;           // ld of X when X=Asub
		size_t ldxc = nr;      // ld of X when X=Csub
		size_t ldy ;           // ld of Y
		size_t la, ca, lb, cb; // lines and columns in A,B sub matrices

		// Three temporary submatrices are required
		typename Field::Element* X = new typename Field::Element[mr*std::max(nr,kr)];
		typename Field::Element* Y = new typename Field::Element[nr*std::max(kr,mr)];

		if (ta == FflasTrans) {
			A21 = A + mr;
			A12 = A + kr*lda;
			A22 = A12 + mr;
			ldxa = mr ;
			la = kr ;
			ca = mr ;
		}
		else { // ta == FflasNoTrans
			A12 = A + kr;
			A21 = A + mr*lda;
			A22 = A21 + kr;
			ldxa = kr ;
			la = mr ;
			ca = kr ;

		}
		if (tb == FflasTrans) {
			B21 = B + kr;
			B12 = B + nr*ldb;
			B22 = B12 + kr;
			ldy = kr ;
			lb = nr ;
			cb = kr ;

		}
		else { // ta == FflasNoTrans
			B12 = B + nr;
			B21 = B + kr*ldb;
			B22 = B21 + nr;
			ldy = nr ;
			lb = kr ;
			cb = nr ;
		}

		// Z1 = C22 - C12         in C22
		fsubin(F,mr,nr,C12,ldc,C22,ldc);
		// Z3 = C12-C21           in C12
		fsubin(F,mr,nr,C21,ldc,C12,ldc);
		// S1 = A21 + A22         in X
		fadd(F,la,ca,A21,lda,A22,lda,X,ldxa);
		// T1 = B12 - B11         in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,ldy);
		// P5 = a S1 T1 + b Z3    in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, beta, C12, ldc, kmax,w-1,base);
		// S2 = S1 - A11          in X
		fsubin(F,la,ca,A11,lda,X,ldxa);
		// T2 = B22 - T1          in Y
		fsub(F,lb,cb,B22,ldb,Y,ldy,Y,ldy);
		// P6 = a S2 T2 + b C21   in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, beta, C21, ldc, kmax,w-1,base);
		// S4 = A12 - S2          in X
		fsub(F,la,ca,A12,lda,X,ldxa,X,ldxa);
		// W1 = P5 + beta Z1      in C22
		fadd(F,mr,nr,C12,ldc,beta,C22,ldc,C22,ldc);
		// P3 = a S4 B22 + P5     in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, B22, ldb, F.zero, Y, nr, kmax,w-1,base);
		fadd(F,mr,nr,Y,nr,C12,ldc,C12,ldc);
		// P1 = a A11 B11         in X
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X, ldxc, kmax,w-1,base);
		// U2 = P6 + P1           in C21
		faddin(F,mr,nr,X,ldxc,C21,ldc);
		// P2 = a A12 B21 + b C11 in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, Y, nr, kmax,w-1,base);
		fadd(F,mr,nr,Y,nr,beta,C11,ldc,C11,ldc);
		// U1 = P1 + P2           in C11
		faddin(F,mr,nr,X,ldxc,C11,ldc);
		// U5 = U2 + P3           in C12
		faddin(F,mr,nr,C21,ldc,C12,ldc);
		// S3 =  A11 - A21        in X ;
		fsub(F,la,ca,A11,lda,A21,lda,X,ldxa);
		// T3 = B22 - B12         in Y
		fsub(F,lb,cb,B22,ldb,B12,ldb,Y,ldy);
		// U3 = a S3 T3 + U2      in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, X, ldxa, Y, ldy, F.one, C21, ldc, kmax,w-1,base);
		// U7 = U3 + W1           in C22
		faddin(F,mr,nr,C21,ldc,C22,ldc);
		// T1_ = B12 - B11        in Y
		fsub(F,lb,cb,B12,ldb,B11,ldb,Y,ldy);
		// T2_ = B22 - T1_        in Y
		fsub(F,lb,cb,B22,ldb,Y,ldy,Y,ldy);
		// T4 = T2_ - B21         in Y
		fsub(F,lb,cb,Y,ldy,B21,ldb,Y,ldy);
		// U6 = -a A22 T4 + U3    in C21;
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, Y, ldy, F.zero, X, nr, kmax,w-1,base);
		delete[] Y;
		fsub(F,mr,nr,C21,ldc,X,nr,C21,ldc);
		delete[] X;


	} // WinogradAccOld


} // BLAS3

} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_acc2_INL

