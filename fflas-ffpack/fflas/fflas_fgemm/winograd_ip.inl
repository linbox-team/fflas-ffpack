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

/** @file fflas/fflas_fgemm/winograd_ip.inl
 * @ingroup MMalgos
 * @brief Winograd implementation
 * @bib ISSAC09 Scheduling
 */

#ifndef __FFLASFFPACK_fgemm_winograd_ip_INL
#define __FFLASFFPACK_fgemm_winograd_ip_INL

namespace FFLAS { namespace BLAS3 {

	template < class Field >
	inline void WinogradIP (const Field& F,
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element alpha,
			      typename Field::Element* A,const size_t lda,
			      typename Field::Element* B,const size_t ldb,
			      const typename Field::Element  beta,
			      typename Field::Element * C, const size_t ldc,
			      const size_t kmax, const size_t w, const FFLAS_BASE base)
	{

		FFLASFFPACK_check(F.isZero(beta));

		// FFLASFFPACK_check(mr == nr && mr == kr);
		FFLASFFPACK_check(kr == nr);

		size_t lb, cb, la, ca;
		size_t ldxa, ldxb;
		typename Field::Element * A11=A, *A12, *A21, *A22;
		typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C21+nr;


		if (ta == FflasTrans) {
			A21  = A + mr;
			A12  = A + kr*lda;
			A22  = A12 + mr;
			la   = kr;
			ca   = mr;
			ldxa = lda;
		}
		else {
			A12  = A + kr;
			A21  = A + mr*lda;
			A22  = A21 + kr;
			la   = mr;
			ca   = kr;
			ldxa = lda;
		}
		if (tb == FflasTrans) {
			B21  = B + kr;
			B12  = B + nr*ldb;
			B22  = B12 + kr;
			lb   = nr;
			cb   = kr;
			ldxb = ldb;
		}
		else {
			B12  = B + nr;
			B21  = B + kr*ldb;
			B22  = B21 + nr;
			lb   = kr;
			cb   = nr;
			ldxb = ldb;
		}


		// S3 = A11 - A21         in C11
		fsub(F,la,ca,A11,ldxa,A21,ldxa,C11,ldc);
		// S1 =  A21 + A22        in A21
		faddin(F,la,ca,A22,ldxa,A21,ldxa);
		// T1 = B12 - B11         in C22
		fsub(F,lb,cb,B12,ldxb,B11,ldxb,C22,ldc);
		// T3 = B22 - B12         in B12
		fsub(F,lb,cb,B22,ldxb,B12,ldxb,B12,ldxb);
		// P7 = S3 T3             in C21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, C11, ldc, B12, ldxb, F.zero, C21, ldc, kmax, w-1, base);
		// S2 = S1 - A11          in C12
		fsub(F,la,ca,A21,ldxa,A11,ldxa,C12,ldc);
		// P1 = A11 B11           in C11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, ldxa, B11, ldxb, F.zero, C11, ldc, kmax, w-1, base);
		// T2 = B22 - T1          in B11
		fsub(F,lb,cb,B22,ldxb,C22,ldc,B11,ldxb);
		// P5 = S1 T1             in A11
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A21, ldxa, C22, ldc, F.zero, A11, ldxa, kmax, w-1, base);
		// T4 = T2 - B21          in C22
		fsub(F,lb,cb,B11,ldxb,B21,ldxb,C22,ldc);
		// P4 = A22 T4            in A21
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, ldxa, C22, ldc, F.zero, A21, ldxa, kmax, w-1, base);
		// S4 = A12 - S2          in A22
		fsub(F,la,ca,A12,ldxa,C12,ldc,A22,ldxa);
		// P6 = S2 T2             in C22
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, C12, ldc, B11, ldxb, F.zero, C22, ldc, kmax, w-1, base);
		// U2 = P1 + P6           in C22
		faddin(F,mr,nr,C11,ldc,C22,ldc);
		// P2 = A12 B21           in C12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, ldxa, B21, ldxb, F.zero, C12, ldc, kmax, w-1, base);
		// U1 = P1 + P2           in C11
		faddin(F,mr,nr,C12,ldc,C11,ldc);
		// U4 = U2 + P5           in C12
		fadd(F,mr,nr,C22,ldc,A11,ldxa,C12,ldc);
		// U3 = U2 + P7           in C22
		faddin(F,mr,nr,C21,ldc,C22,ldc);
		// U6 = U3 - P4           in C21
		fsub(F,mr,nr,C22,ldc,A21,ldxa,C21,ldc);
		// U7 = U3 + P5           in C22
		faddin(F,mr,nr,A11,ldxa,C22,ldc);
		// P3 = S4 B22            in A12
		Protected::WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, ldxa, B22, ldxb, F.zero, A12, ldxa, kmax, w-1, base);
		// U5 = U4 + P3           in C12
		faddin(F,mr,nr,A12,ldxa,C12,ldc);



	} // WinogradIP

} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_winograd_ip_INL

