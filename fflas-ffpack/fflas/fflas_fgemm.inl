/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fgemm.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet  < Clement.Pernet@imag.fr >
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

#ifndef __FFLASFFPACK_fgemm_INL
#define __FFLASFFPACK_fgemm_INL


namespace FFLAS {

	namespace Protected {

		// Note:
		// The domain is supposed to be a field since some divisions are required for efficiency purposes
		// An alternative has to be written for finite rings if necessary

		// Classic Multiplication over double
		// Classic multiplication over a finite field
		template  < class Field >
		inline void ClassicMatmul (const Field& F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const typename Field::Element alpha,
					   const typename Field::Element * A, const size_t lda,
					   const typename Field::Element * B, const size_t ldb,
					   const typename Field::Element beta,
					   typename Field::Element* C, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			typename Field::Element tmp;

			size_t k2 = std::min(k,kmax); // Size of the blocks

			if (k2 > 1) {
				if (base == FflasDouble){
					DoubleDomain::Element alphad, betad;
					DoubleDomain::Element * Add = new DoubleDomain::Element[m*k2];
					DoubleDomain::Element * Bdd = new DoubleDomain::Element[k2*n];
					DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];

					size_t nblock = k / kmax;
					size_t remblock = k % kmax;
					if (!remblock) {
						remblock = kmax ;
						--nblock;
					}
					if (F.isMOne( beta)) betad = -1.0;
					else F.convert (betad, beta);

					if (F.isMOne( alpha)) alphad = -1.0;
					else {
						alphad = 1.0;
						if (! F.isOne( alpha)) {
							// Compute y = A*x + beta/alpha.y
							// and after y *= alpha
							FFLASFFPACK_check(!F.isZero(alpha));
							F.div (tmp, beta, alpha);
							F.convert (betad, tmp);
						}
					}

					size_t dlda, dldb;
					if (!F.isZero(beta))
						MatF2MatD (F, Cd, n, C, ldc, m, n);

					if (ta == FflasTrans) {
						dlda = m;
						MatF2MatD (F, Add, dlda, A+k2*nblock*lda, lda, remblock, m);
					}
				       	else {
						dlda = k2;
						MatF2MatD (F, Add, dlda, A+k2*nblock, lda, m, remblock);
					}
					if (tb == FflasTrans) {
						dldb = k2;
						MatF2MatD (F, Bdd, k2, B+k2*nblock, ldb, n, remblock);
					}
				       	else {
						dldb = n;
						MatF2MatD (F, Bdd, dldb, B+k2*nblock*ldb, ldb, remblock, n);
					}

					ClassicMatmul (DoubleDomain(), ta, tb, m, n, remblock, alphad, Add, dlda,
						       Bdd, dldb, betad, Cd, n, kmax,base );

					MatD2MatF (F, C, ldc, Cd, n, m, n);
					MatF2MatD (F, Cd, n, C, ldc, m, n);

					for (size_t i = 0; i < nblock; ++i) {
						if (ta == FflasTrans) MatF2MatD (F, Add, dlda, A+k2*i*lda, lda, k2, m);
						else MatF2MatD (F, Add, dlda,  A+k2*i, lda, m, k2);

						if (tb == FflasTrans) MatF2MatD (F, Bdd, dldb, B+k2*i, ldb, n, k2);
						else MatF2MatD (F, Bdd, dldb, B+k2*i*ldb, ldb, k2, n);

						ClassicMatmul (DoubleDomain(), ta, tb, m, n, k2, alphad, Add, dlda,
							       Bdd, dldb, 1.0, Cd, n, kmax,base);
						MatD2MatF (F, C, ldc, Cd, n, m, n);
						MatF2MatD (F, Cd, n, C, ldc, m, n);
					}
					if ((!F.isOne( alpha)) && (!F.isMOne( alpha)))
						fscalin(F,m,n,alpha,C,ldc);
					delete[] Add;
					delete[] Bdd;
					delete[] Cd;
				}
				else {
					FloatDomain::Element alphad, betad;
					FloatDomain::Element * Add = new FloatDomain::Element[m*k2];
					FloatDomain::Element * Bdd = new FloatDomain::Element[k2*n];
					FloatDomain::Element * Cd = new FloatDomain::Element[m*n];

					size_t nblock = k / kmax;
					size_t remblock = k % kmax;
					if (!remblock) {
						remblock = kmax;
						--nblock;
					}
					if (F.isMOne( beta)) betad = -1.0;
					else F.convert (betad, beta);

					if (F.isMOne( alpha)) alphad = -1.0;
					else {
						alphad = 1.0;
						if (! F.isOne( alpha)) {
							// Compute y = A*x + beta/alpha.y
							// and after y *= alpha
							FFLASFFPACK_check(!F.isZero(alpha));
							F.div (tmp, beta, alpha);
							F.convert (betad, tmp);
						}
					}

					size_t dlda, dldb;
					if (!F.isZero(beta))
						MatF2MatFl (F, Cd, n, C, ldc, m, n);

					if (ta == FflasTrans) {
						dlda = m;
						MatF2MatFl (F, Add, dlda, A+k2*nblock*lda, lda, remblock, m);
					}
				       	else {
						dlda = k2;
						MatF2MatFl (F, Add, dlda, A+k2*nblock, lda, m, remblock);
					}
					if (tb == FflasTrans) {
						dldb = k2;
						MatF2MatFl (F, Bdd, k2, B+k2*nblock, ldb, n, remblock);
					}
				       	else {
						dldb = n;
						MatF2MatFl (F, Bdd, dldb, B+k2*nblock*ldb, ldb, remblock, n);
					}

					ClassicMatmul (FloatDomain(), ta, tb, m, n, remblock, alphad, Add, dlda,
						       Bdd, dldb, betad, Cd, n, kmax,base );
					MatFl2MatF (F, C, ldc, Cd, n, m, n);
					MatF2MatFl (F, Cd, n, C, ldc, m, n);
					for (size_t i = 0; i < nblock; ++i) {
						if (ta == FflasTrans) MatF2MatFl (F, Add, dlda, A+k2*i*lda, lda, k2, m);
						else MatF2MatFl (F, Add, dlda,  A+k2*i, lda, m, k2);
						if (tb == FflasTrans) MatF2MatFl (F, Bdd, dldb, B+k2*i, ldb, n, k2);
						else MatF2MatFl (F, Bdd, dldb, B+k2*i*ldb, ldb, k2, n);

						ClassicMatmul (FloatDomain(), ta, tb, m, n, k2, alphad, Add, dlda,
							       Bdd, dldb, 1.0, Cd, n, kmax,base);
						MatFl2MatF (F, C, ldc, Cd, n, m, n);
						MatF2MatFl (F, Cd, n, C, ldc, m, n);
					}
					if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
						fscalin(F,m,n,alpha,C,ldc);
					}
					delete[] Add;
					delete[] Bdd;
					delete[] Cd;
				}
			}
			else { // k2 == 1
				// Standard algorithm is performed over the Field, without conversion
				if (F.isZero (beta))
					for (size_t i = 0; i < m; ++i)
						for (size_t j = 0; j < n; ++j)
							F.assign (*(C+i*ldc+j), F.zero);
				else {
					typename Field::Element betadivalpha;
					FFLASFFPACK_check(!F.isZero(alpha));
					F.div (betadivalpha, beta, alpha);
					fscalin(F,m,n,betadivalpha,C,ldc);
				}
				if (ta == FflasNoTrans)
					if (tb == FflasNoTrans)
						for (size_t i = 0; i < m; ++i)
							for (size_t l = 0; l < k; ++l)
								for (size_t j = 0; j < n; ++j)
									F.axpyin (*(C+i*ldc+j), *(A+i*lda+l), *(B+l*ldb+j));
					else
						for (size_t i = 0; i < m; ++i)
							for (size_t j = 0; j < n; ++j)
								for (size_t l = 0; l < k; ++l)
									F.axpyin (*(C+i*ldc+j), *(A+i*lda+l), *(B+j*ldb+l));
				else
					if (tb == FflasNoTrans)
						for (size_t i = 0; i < m; ++i)
							for (size_t l = 0; l < k; ++l)
								for (size_t j = 0; j < n; ++j)
									F.axpyin (*(C+i*ldc+j), *(A+l*lda+i), *(B+l*ldb+j));
					else
						for (size_t i = 0; i < m; ++i)
							for (size_t j = 0; j < n; ++j)
								for (size_t l = 0; l < k; ++l)
									F.axpyin (*(C+i*ldc+j), *(A+l*lda+i), *(B+j*ldb+l));
				if (! F.isOne(alpha))
					fscalin(F,m,n,alpha,C,ldc);
			}
		}

		template<>
		inline void ClassicMatmul (const DoubleDomain& ,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const DoubleDomain::Element alpha,
					   const DoubleDomain::Element * Ad, const size_t lda,
					   const DoubleDomain::Element * Bd, const size_t ldb,
					   const DoubleDomain::Element beta,
					   DoubleDomain::Element * Cd, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{

			cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
				     (int)m, (int)n, (int)k, (DoubleDomain::Element) alpha,
				     Ad, (int)lda, Bd, (int)ldb, (DoubleDomain::Element) beta,Cd, (int)ldc);
		}

		template  <>
		inline void ClassicMatmul (const FloatDomain& F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const FloatDomain::Element alpha,
					   const FloatDomain::Element * Ad, const size_t lda,
					   const FloatDomain::Element * Bd, const size_t ldb,
					   const FloatDomain::Element beta,
					   FloatDomain::Element * Cd, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			cblas_sgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
				     (int)m, (int)n, (int)k, (FloatDomain::Element) alpha,
				     Ad, (int)lda, Bd, (int)ldb, (FloatDomain::Element) beta,Cd, (int)ldc);
		}

		template <>
		inline void ClassicMatmul (const FFPACK:: ModularBalanced<double> & F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const double alpha,
					   const double * A, const size_t lda,
					   const double * B, const size_t ldb,
					   const double beta,
					   double* C, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			double _alpha, _beta;
			// To ensure the initial computation with beta
			size_t k2 = std::min(k,kmax);
			size_t nblock = k / kmax;
			size_t remblock = k % kmax;
			if (!remblock) {
				remblock = kmax;
				--nblock;
			}
			if (F.isMOne( beta)) _beta = -1.0;
			else _beta = beta;
			if (F.isMOne( alpha)) _alpha = -1.0;
			else{
				_alpha = 1.0;
				if (! F.isOne( alpha)) {
					// Compute y = A*x + beta/alpha.y
					// and after y *= alpha
					FFLASFFPACK_check(!F.isZero(alpha));
					F.divin (_beta, alpha);
				}
			}
			size_t shiftA, shiftB;
			if (ta == FflasTrans) shiftA = k2*lda;
			else shiftA = k2;
			if (tb == FflasTrans) shiftB = k2;
			else shiftB = k2*ldb;

			ClassicMatmul (DoubleDomain(), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
				       B+nblock*shiftB, ldb, _beta, C, ldc, kmax,base );
			finit(F,m,n,C,ldc);
			for (size_t i = 0; i < nblock; ++i) {
				ClassicMatmul (DoubleDomain(), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
					       B+i*shiftB, ldb, F.one, C, ldc, kmax,base);
				finit(F,m,n,C,ldc);
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
		}


		template <>
		inline void ClassicMatmul (const  FFPACK:: ModularBalanced<float> & F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const float alpha,
					   const float * A, const size_t lda,
					   const float * B, const size_t ldb,
					   const float beta,
					   float* C, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			float _alpha, _beta;
			// To ensure the initial computation with beta
			size_t k2 = std::min(k,kmax);
			size_t nblock = k / kmax;
			size_t remblock = k % kmax;
			if (!remblock) {
				remblock = kmax;
				--nblock;
			}
			if (F.isMOne( beta)) _beta = -1.0;
			else _beta = beta;
			if (F.isMOne( alpha)) _alpha = -1.0;
			else{
				_alpha = 1.0;
				if (! F.isOne( alpha)) {
					// Compute y = A*x + beta/alpha.y
					// and after y *= alpha
					FFLASFFPACK_check(!F.isZero(alpha));
					F.divin (_beta, alpha);
				}
			}
			size_t shiftA, shiftB;
			if (ta == FflasTrans) shiftA = k2*lda;
			else shiftA = k2;
			if (tb == FflasTrans) shiftB = k2;
			else shiftB = k2*ldb;

			ClassicMatmul (FloatDomain(), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
				       B+nblock*shiftB, ldb, _beta, C, ldc, kmax,base);
			finit(F,m,n,C,ldc);
			for (size_t i = 0; i < nblock; ++i) {
				ClassicMatmul (FloatDomain(), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
					       B+i*shiftB, ldb, F.one, C, ldc, kmax,base);
				finit(F,m,n,C,ldc);
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
		}


		template <>
		inline void ClassicMatmul (const  FFPACK:: Modular<double> & F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const double alpha,
					   const double * A, const size_t lda,
					   const double * B, const size_t ldb,
					   const double beta,
					   double* C, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			double _alpha, _beta;
			// To ensure the initial computation with beta
			size_t k2 = std::min(k,kmax);
			size_t nblock = k / kmax;
			size_t remblock = k % kmax;
			if (!remblock) {
				remblock = kmax;
				--nblock;
			}
			if (F.isMOne( beta))
				_beta = -1.0;
			else
				_beta = beta;
			if (F.isMOne( alpha))
				_alpha = -1.0;
			else {
				_alpha = 1.0;
				if (! F.isOne( alpha)) {
					// Compute y = A*x + beta/alpha.y
					// and after y *= alpha
					FFLASFFPACK_check(!F.isZero(alpha));
					F.divin (_beta, alpha);
				}
			}
			size_t shiftA, shiftB;
			if (ta == FflasTrans)
				shiftA = k2*lda;
			else
				shiftA = k2;
			if (tb == FflasTrans)
				shiftB = k2;
			else
				shiftB = k2*ldb;

			ClassicMatmul (DoubleDomain(), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
				       B+nblock*shiftB, ldb, _beta, C, ldc, kmax, base );

#if 0 /* timing */
			Timer t;
			t.clear();
			t.start();
#endif
			//#pragma omp parallel for schedule(static) private (Ci)
			//!@todo init only if remblock!=0 and _beta == 0
			finit(F,m,n,C,ldc);
#if 0
			t.stop();
			std::cerr<<"Reduction dans Classic -> "<<t.realtime()<<std::endl;
#endif


			for (size_t i = 0; i < nblock; ++i) {
				ClassicMatmul (DoubleDomain(), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
					       B+i*shiftB, ldb, F.one, C, ldc, kmax, base);
				//#pragma omp parallel for schedule(static) private (Ci)
				finit(F,m,n,C,ldc);
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
		}

		template <>
		inline void ClassicMatmul (const  FFPACK:: Modular<float> & F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m, const size_t n,const size_t k,
					   const float alpha,
					   const float * A, const size_t lda,
					   const float * B, const size_t ldb,
					   const float beta,
					   float* C, const size_t ldc,
					   const size_t kmax, const FFLAS_BASE base)
		{
			float _alpha, _beta;
			// To ensure the initial computation with beta
			size_t k2 = std::min(k,kmax);
			size_t nblock = k / kmax;
			size_t remblock = k % kmax;
			if (!remblock) {
				remblock = kmax;
				--nblock;
			}
			if (F.isMOne( beta)) _beta = -1.0;
			else _beta = beta;
			if (F.isMOne( alpha)) _alpha = -1.0;
			else{
				_alpha = 1.0;
				if (! F.isOne( alpha)) {
					// Compute y = A*x + beta/alpha.y
					// and after y *= alpha
					FFLASFFPACK_check(!F.isZero(alpha));
					F.divin (_beta, alpha);
				}
			}
			size_t shiftA, shiftB;
			if (ta == FflasTrans) shiftA = k2*lda;
			else shiftA = k2;
			if (tb == FflasTrans) shiftB = k2;
			else shiftB = k2*ldb;

			ClassicMatmul (FloatDomain(), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
				       B+nblock*shiftB, ldb, _beta, C, ldc, kmax,base);
			finit(F,m,n,C,ldc);
			for (size_t i = 0; i < nblock; ++i) {
				ClassicMatmul (FloatDomain(), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
					       B+i*shiftB, ldb, F.one, C, ldc, kmax,base);
				finit(F,m,n,C,ldc);
				// for (float * Ci = C; Ci != C+m*ldc; Ci += ldc)
					// for (size_t j=0; j < n;++j)
						// F.init(*(Ci+j),*(Ci+j));
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
		}


	// Winograd Multiplication  A(n*k) * B(k*m) in C(n*m)
	// Computation of the 22 Winograd's operations
	template < class Field >
	inline void WinoCalc (const Field& F,
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element alpha,
			      const typename Field::Element* A,const size_t lda,
			      const typename Field::Element* B,const size_t ldb,
			      const typename Field::Element beta,
			      typename Field::Element * C, const size_t ldc,
			      const size_t kmax, const size_t w, const FFLAS_BASE base)
	{

		typename Field::Element mbeta;
		F.neg(mbeta,beta);
		size_t imaxb, jmaxb, imaxa, jmaxa, ldx2;
		size_t x3rd = std::max(mr,kr);
		const typename Field::Element* d11,*d12,*d21,*d22;
		typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2,*dx3;
		const typename Field::Element * A11=A, *A12, *A21, *A22;
		const typename Field::Element * B11=B, *B12, *B21, *B22;
		typename Field::Element * C11=C, *C12=C+nr, *C21=C+mr*ldc, *C22=C+nr+mr*ldc;


		if (F.isZero(beta)){
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
			typename Field::Element* X1 = new typename Field::Element[mr*x1rd];		// S3 = A11 - A21 in X1
			d11 = A11; d21 = A21; dx1 = X1;
			fsub(F,imaxa,jmaxa,d11,lda,d21,lda,dx1,ldx1);

			// P7 = alpha . S3 * T3  in C21
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C21, ldc, kmax, w-1, base);

			// T1 = B12 - B11 in X2
			d11 = B11; d12 = B12; dx2 = X2;
			fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx2,ldx2);

			// S1 = A21 + A22 in X1

			d21 = A21; d22 = A22; dx1 = X1;
			fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx1,ldx1);

			// P5 = alpha . S1*T1 in C22
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C22, ldc, kmax, w-1, base);

			// T2 = B22 - T1 in X2
			d22 = B22; dx2 = X2;
			fsub(F,imaxb,jmaxb,d22,ldb,dx2,ldx2,dx2,ldx2);

			// S2 = S1 - A11 in X1
			d11 = A11; dx1 = X1;
			fsubin(F,imaxa,jmaxa,d11,lda,dx1,ldx1);

			// P6 = alpha . S2 * T2 in C12
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, X2, ldx2, F.zero, C12, ldc, kmax, w-1, base);

			// S4 = A12 -S2 in X1
			d12 = A12; dx1 = X1;
			fsub(F,imaxa,jmaxa,d12,lda,dx1,ldx1,dx1,ldx1);

			// P3 = alpha . S4*B22 in C11
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X1, ldx1, B22, ldb, F.zero, C11, ldc, kmax, w-1, base);

			// P1 = alpha . A11 * B11 in X1
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, kmax, w-1, base);



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
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X2, ldx2, F.zero, C11, ldc, kmax, w-1, base);

			delete[] X2;
			// U6 = U3 - P4 in C21
			d21c = C21; d11c = C11;
			fsubin(F,mr,nr,d11c,ldc,d21c,ldc);

			// P2 = alpha . A12 * B21  in C11
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, F.zero, C11, ldc, kmax,w-1, base);

			//  U1 = P2 + P1 in C11
			d11c = C11; dx1 = X1;
			faddin(F,mr,nr,dx1,nr,d11c,ldc);

			delete[] X1;

		}
	       	else { // beta non zero
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

#ifdef NEWWINO
#if 0
			std::cerr<<"New Wino"<<std::endl;
			// C22 = C22 - C12
			d12c = C12;
			d22c = C22;
			for (size_t i = 0; i <  mr; ++i, d12c += ldc, d22c += ldc)
				for (size_t j = 0; j < nr; ++j)
					F.subin (*(d22c + j), *(d12c + j));
#endif


			// T1 = B12 - B11 in X3
			d11 = B11; d12 = B12; dx3 = X3;
			fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx3,ldx3);

			// S1 = A21 + A22 in X2
			d21 = A21; d22 = A22; dx2 = X2;
			fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx2,ldx2);

			// P5 = alpha . S1*T1 + beta . C12 in C12
			//WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, beta, C12, ldc, kmax, w-1,base);
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.zero, X1, nr, kmax, w-1,base);

			// C22 = P5 + beta C22 in C22
			d22c = C22; dx1 = X1;
			for (size_t i = 0; i < mr; ++i, dx1 += nr, d22c += ldc)
				for (size_t j=0;j < nr;++j) {
					//! @todo can merge ops ?
					F.mulin (*(d22c + j), beta);
					F.addin (*(d22c + j), *(dx1 + j));
				}

			// C12 = P5 + beta C12 in C12
			dx1 = X1; d12c = C12;
			for (size_t i = 0; i < mr; ++i, d12c += ldc, dx1 += nr)
				for (size_t j=0;j < nr;++j) {
					//! @todo can merge ops ?
					F.mulin (*(d12c + j), beta);
					F.addin (*(d12c + j), *(dx1 + j));
				}

			// P1 = alpha . A11 * B11 in X1
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X1, nr, kmax, w-1,base);


			// P2 = alpha . A12 * B21 + beta . C11  in C11
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, kmax,w-1,base);

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
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.one, X1, nr, kmax, w-1,base);

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
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldx3, mbeta, C21, ldc, kmax, w-1,base);

			// U5 = P3 + U4 = alpha . S4*B22 + U4 in C12
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, B22, ldb, F.one, C12, ldc, kmax, w-1,base);

			// T3 = B22 - B12 in X3
			d12 = B12; d22 = B22; dx3 = X3;
			fsub(F,imaxb,jmaxb,d22,ldb,d12,ldb,dx3,ldx3);

			// S3 = A11 - A21 in X2
			d11 = A11; d21 = A21; dx2 = X2;
			fsub(F,imaxa,jmaxa,d11,lda,d21,lda,dx2,ldx2);

			// U3 = P7 + U2  = alpha . S3 * T3 + U2 in X1
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.one, X1, nr, kmax, w-1,base);

			// U7 =  U3 + C22 in C22
			d22c = C22; dx1 = X1; d12c = C12;
			faddin(F,mr,nr,dx1,nr,d22c,ldc);

			// U6 = U3 - P4 in C21
			dx1 = X1; d21c = C21;
			fsub(F,mr,nr,dx1,nr,d21c,ldc,d21c,ldc);
#else
			// P2 = alpha . A12 * B21 + beta . C11  in C11
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, beta, C11, ldc, kmax,w-1,base);

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
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, beta, C22, ldc, kmax, w-1,base);

			// T1 = B12 - B11 in X3
			d11 = B11; d12 = B12; dx3 = X3;
			fsub(F,imaxb,jmaxb,d12,ldb,d11,ldb,dx3,ldx3);

			// S1 = A21 + A22 in X2
			d21 = A21; d22 = A22; dx2 = X2;
			fadd(F,imaxa,jmaxa,d21,lda,d22,lda,dx2,ldx2);

			// P5 = alpha . S1*T1 + beta . C12 in C12
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, beta, C12, ldc, kmax, w-1,base);

			// T2 = B22 - T1 in X3
			d22 = B22; dx3 = X3;
			fsub(F,imaxb,jmaxb,d22,ldb,dx3,ldx3,dx3,ldx3);

			// S2 = S1 - A11 in X2
			d11 = A11; dx2 = X2;
			fsubin(F,imaxa,jmaxa,d11,lda,dx2,ldx2);

			// P6 = alpha . S2 * T2 in X1
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, X3, ldx3, F.zero, X1, nr, kmax, w-1,base);

			// T4 = T2 - B21 in X3
			d21 = B21;dx3=X3;
			fsubin(F,imaxb,jmaxb,d21,ldb,dx3,ldx3);

			// S4 = A12 -S2 in X2
			d12 = A12; dx2 = X2;
			fsub(F,imaxa,jmaxa,d12,lda,dx2,ldx2,dx2,ldx2);

			// P4 = alpha . A22 * T4 - beta . C21 in C21
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A22, lda, X3, ldx3, mbeta, C21, ldc, kmax, w-1,base);

			// P1 = alpha . A11 * B11 in X3
			WinoMain (F, ta, tb, mr, nr, kr, alpha, A11, lda, B11, ldb, F.zero, X3, nr, kmax, w-1,base);

			//  U1 = P2 + P1 in C11
			d11c = C11; dx3 = X3;
			faddin(F,mr,nr,dx3,ldx3,d11c,ldc);

			// U2 = P1 + P6 in tmpU2  and
			// U3 = P7 + U2 in tmpU3  and
			// U7 = P5 + U3 in C22    and
			// U4 = P5 + U2 in C12    and
			// U6 = U3 - P4 in C21    and
			typename Field::Element tmpU2, tmpU3;
			d12c = C12; dx1=X1; dx3=X3; d21c = C21; d22c = C22;
			for (size_t i = 0; i < mr;
			     ++i, d12c += ldc, dx1 += nr, dx3 += nr, d22c+=ldc, d21c += ldc) {
				for (size_t j=0;j < nr;++j) {
					F.add (tmpU2, *(dx3 + j), *(dx1 + j));    // temporary U2 = P1 + P6
					F.add (tmpU3, tmpU2, *(d22c + j));      // temporary U3 = U2 + P7
					F.add (*(d22c + j), *(d12c + j), tmpU3);  // U7 = P5 + U3 in C22
					F.addin (*(d12c + j), tmpU2);             // U4 = P5 + U2 in C12
					F.sub (*(d21c + j), tmpU3, *(d21c + j)); // U6 = U3 - P4 in C21
				}
			}

			delete[] X1;
			delete[] X3;

			// P3 = alpha . S4*B22 in X1
			WinoMain (F, ta, tb, mr, nr, kr, alpha, X2, ldx2, B22, ldb, F.one, C12, ldc, kmax, w-1,base);

			// U5 = P3 + U4 in C12
#if 0
			d12c = C12; dx1 = X1;
			for (size_t i = 0; i < mr; ++i, d12c += ldc, dx1 += nr)
				for (size_t j = 0; j < nr; ++j)
					F.addin (*(d12c + j), *(dx1 + j));
#endif
#endif
			delete[] X2;
		}
	}


		// Control the switch with classic multiplication
		// Fix-up for odd-sized matrices using dynamic pealing
		// for matrices over double
		template <>
		inline  void WinoMain (const DoubleDomain& D,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n, const size_t k,
				       const DoubleDomain::Element alpha,
				       const DoubleDomain::Element * A, const size_t lda,
				       const DoubleDomain::Element * B, const size_t ldb,
				       const DoubleDomain::Element beta,
				       DoubleDomain::Element * C, const size_t ldc,
				       const size_t kmax, const size_t w, const FFLAS_BASE base)
		{

			if (w == 0)
				ClassicMatmul (D, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			else{
				WinoCalc (D, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb,
					  beta, C, ldc, kmax, w,base);
				DynamicPealing (D, ta, tb, m, n, k, alpha, A, lda, B, ldb,
						beta, C, ldc, kmax);
			}
		}

		template <>
		inline  void WinoMain (const FloatDomain& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n, const size_t k,
				       const FloatDomain::Element alpha,
				       const FloatDomain::Element * A, const size_t lda,
				       const FloatDomain::Element * B, const size_t ldb,
				       const FloatDomain::Element beta,
				       FloatDomain::Element * C, const size_t ldc,
				       const size_t kmax, const size_t w, const FFLAS_BASE base)
		{

			if (w == 0) {
				ClassicMatmul (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			}
			else{
				WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb,
					  beta, C, ldc, kmax, w,base);
				DynamicPealing (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
						beta, C, ldc, kmax);
			}
		}

		template <class Field>
		inline void  WinoMain (const Field& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n, const size_t k,
				       const typename Field::Element alpha,
				       const typename Field::Element* A,const size_t lda,
				       const typename Field::Element* B,const size_t ldb,
				       const typename Field::Element beta,
				       typename Field::Element * C, const size_t ldc,
				       const size_t kmax, const size_t w,
				       const FFLAS_BASE base)
		{


			if (w == 0) // Winograd - >  Classic
				ClassicMatmul (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			else {
				if (k <= kmax) { // switch on floating point
					if (base == FflasDouble){
						DoubleDomain::Element alphad, betad;
						typename Field::Element _betabis;

						if (F.isMOne(alpha)) {
							alphad = -1.0;
							F.convert (betad, beta);
						}
					       	else {
							if (! F.isOne ( alpha)) {
								// Compute C = A*B + beta/alpha.C
								// and after C *= alpha
								FFLASFFPACK_check(!F.isZero(alpha));
								F.div (_betabis, beta, alpha);
								F.convert (betad, _betabis);
							}
							else
								F.convert (betad, beta);
							alphad = 1.0;
						}
						DoubleDomain::Element * Ad = new DoubleDomain::Element[m*k];
						DoubleDomain::Element * Bd = new DoubleDomain::Element[k*n];
						DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
						// Conversion GFq = >  double
						size_t ma, ka, kb, nb; //mb, na
						if (ta == FflasTrans) { ma = k; ka = m; }
						else { ma = m; ka = k; }
						if (tb == FflasTrans) { kb = n; nb = k; }
						else {  kb = k; nb = n; }

						MatF2MatD (F, Ad, ka, A, lda, ma, ka);
						MatF2MatD (F, Bd, nb, B, ldb, kb, nb);
						if (!F.isZero(beta))
							MatF2MatD (F, Cd, n, C, ldc, m, n);
						// recursive call
						WinoMain (DoubleDomain(), ta, tb, m, n, k, alphad,
							  Ad, ka, Bd, nb, betad, Cd, n, kmax, w,base);
						// Conversion double = >  GFq
						MatD2MatF (F, C, ldc, Cd, n, m, n);

						if (!F.isOne(alpha) &&
						    !F.isMOne (alpha)) {
							// Fix-up: compute C *= alpha
							fscalin(F,m,n,alpha,C,ldc);
						}
						// Temporary double matrices destruction
						delete[] Ad;
						delete[] Bd;
						delete[] Cd;
					}
				       	else {
						FloatDomain::Element alphad, betad;
						typename Field::Element _betabis;

						if (F.isMOne (alpha)) {
							alphad = -1.0;
							F.convert (betad, beta);
						}
					       	else {
							if (! F.isOne ( alpha)) {
								// Compute C = A*B + beta/alpha.C
								// and after C *= alpha
								FFLASFFPACK_check(!F.isZero(alpha));
								F.div (_betabis, beta, alpha);
								F.convert (betad, _betabis);
							}
							else
								F.convert (betad, beta);
							alphad = 1.0;
						}
						FloatDomain::Element * Ad = new FloatDomain::Element[m*k];
						FloatDomain::Element * Bd = new FloatDomain::Element[k*n];
						FloatDomain::Element * Cd = new FloatDomain::Element[m*n];
						// Conversion GFq = >  double
						size_t ma, ka, kb, nb; //mb, na
						if (ta == FflasTrans) { ma = k; ka = m; }
						else { ma = m; ka = k; }
						if (tb == FflasTrans) { kb = n; nb = k; }
						else {  kb = k; nb = n; }

						MatF2MatFl (F, Ad, ka, A, lda, ma, ka);
						MatF2MatFl (F, Bd, nb, B, ldb, kb, nb);
						if (!F.isZero(beta))
							MatF2MatFl (F, Cd, n, C, ldc, m, n);
						// recursive call
						WinoMain (FloatDomain(), ta, tb, m, n, k, alphad,
							  Ad, ka, Bd, nb, betad, Cd, n, kmax, w,base);
						// Conversion double = >  GFq
						MatFl2MatF (F, C, ldc, Cd, n, m, n);

						if (!F.isOne (alpha) &&
						    !F.isMOne (alpha)) {
							// Fix-up: compute C *= alpha
							fscalin(F,m,n,alpha,C,ldc);
						}
						// Temporary double matrices destruction
						delete[] Ad;
						delete[] Bd;
						delete[] Cd;
					}
				}
			       	else{
					WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb,
						  beta, C, ldc, kmax,w,base);
					DynamicPealing (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
							beta, C, ldc, kmax);
				}
			}
		}

		template <>
		inline void WinoMain (const FFPACK:: ModularBalanced<double>& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t m, const size_t n, const size_t k,
				      const double alpha,
				      const double* A, const size_t lda,
				      const double* B, const size_t ldb,
				      const double beta,
				      double * C, const size_t ldc,
				      const size_t kmax, const size_t w, const FFLAS_BASE base)
		{
			if (w == 0)
				ClassicMatmul (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			else {
				if (k <= kmax) { // switch on delayed modulus
					DoubleDomain::Element _alpha, _beta;
					_beta = beta;
					if (F.isMOne( alpha)) _alpha = -1.0;
					else{
						// Compute C = A*B + beta/alpha.C
						// and then C *= alpha
						if (! F.isOne( alpha)) {
							FFLASFFPACK_check(!F.isZero(alpha));
							F.divin (_beta, alpha);
						}
						_alpha = 1.0;
					}

#if 0			/* BB : useless */
					size_t ma, ka, kb, nb; //mb, na;
					if (ta == FflasTrans) { ma = k; ka = m; }
					else { ma = m; ka = k; }
					if (tb == FflasTrans) { kb = n; nb = k; }
					else {  kb = k; nb = n; }
#endif
					// recursive call
					WinoMain (DoubleDomain(), ta, tb, m, n, k, _alpha,
						  A, lda, B, ldb, _beta, C, ldc, kmax, w,base);
					// Modular reduction
					finit(F,m,n,C,ldc);
					if (!F.isOne( alpha) &&
					    !F.isMOne( alpha))
						// Fix-up: compute C *= alpha
						fscalin(F,m,n,alpha,C,ldc);
				}
			       	else {
					WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha,
						  A, lda, B, ldb, beta, C, ldc, kmax,w,base);
					DynamicPealing (F, ta, tb, m, n, k, alpha,
							A, lda, B, ldb, beta, C, ldc, kmax);
				}
			}
		}


		template <>
		inline void WinoMain (const FFPACK:: ModularBalanced<float>& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t m, const size_t n, const size_t k,
				      const float alpha,
				      const float* A, const size_t lda,
				      const float* B, const size_t ldb,
				      const float beta,
				      float * C, const size_t ldc,
				      const size_t kmax, const size_t w, const FFLAS_BASE base)
		{
			if (w == 0)
				ClassicMatmul (F, ta, tb, m, n, k, alpha,
					       A, lda, B, ldb, beta, C, ldc, kmax,base);
			else {
				if (k <= kmax) { // switch on float
					// Temporary float matrices
					FloatDomain::Element _alpha, _beta;
					_beta = beta;
					if (F.isMOne( alpha)) _alpha = -1.0;
					else {
						// Compute C = A*B + beta/alpha.C
						// and then C *= alpha
						if (! F.isOne( alpha)) {
							FFLASFFPACK_check(!F.isZero(alpha));
							F.divin (_beta, alpha);
						}
						_alpha = 1.0;
					}
#if 0 /* BB : inutile */
					size_t ma, ka, kb, nb; //mb, na;
					if (ta == FflasTrans) { ma = k; ka = m; }
					else { ma = m; ka = k; }
					if (tb == FflasTrans) { kb = n; nb = k; }
					else {  kb = k; nb = n; }
#endif
					// recursive call
					WinoMain (FloatDomain(), ta, tb, m, n, k, _alpha,
						  A, lda, B, ldb, _beta, C, ldc, kmax, w,base);
					// Conversion float = >  GFq
					finit(F,m,n,C,ldc);
					if (!F.isOne( alpha) &&
					    !F.isMOne( alpha))
						// Fix-up: compute C *= alpha
						fscalin(F,m,n,alpha,C,ldc);
				}
			       	else{
					WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha,
						  A, lda, B, ldb, beta, C, ldc, kmax,w,base);
					DynamicPealing (F, ta, tb, m, n, k, alpha,
							A, lda, B, ldb, beta, C, ldc, kmax);
				}
			}
		}

		template <>
		inline void WinoMain (const FFPACK:: Modular<double>& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t m, const size_t n, const size_t k,
				      const double alpha,
				      const double* A, const size_t lda,
				      const double* B, const size_t ldb,
				      const double beta,
				      double * C, const size_t ldc,
				      const size_t kmax, const size_t w, const FFLAS_BASE base)
		{
			if (w == 0)
				ClassicMatmul (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			else {
				if (k <= kmax) { // switch on delayed modulus
					DoubleDomain::Element _alpha, _beta;
					_beta = beta;
					if (F.isMOne( alpha)) _alpha = -1.0;
					else{
						// Compute C = A*B + beta/alpha.C
						// and then C *= alpha
						if (! F.isOne( alpha)) {
							FFLASFFPACK_check(!F.isZero(alpha));
							F.divin (_beta, alpha);
						}
						_alpha = 1.0;
					}

#if 0 /* BB: inutile */
					size_t ma, ka, kb, nb; //mb, na;
					if (ta == FflasTrans) { ma = k; ka = m; }
					else { ma = m; ka = k; }
					if (tb == FflasTrans) { kb = n; nb = k; }
					else {  kb = k; nb = n; }
#endif
					// recursive call
					WinoMain (DoubleDomain(), ta, tb, m, n, k, _alpha,
						  A, lda, B, ldb, _beta, C, ldc, kmax, w,base);
					// Modular reduction
#if 0 /*  timing */
					Timer t;
					t.clear();
					t.start();
#endif

					// #pragma omp parallel for schedule(static) private (Ci)
					finit(F,m,n,C,ldc);
#if 0
					t.stop();
					std::cerr<<"Reduction -> "<<t.realtime()<<std::endl;
#endif

					if (!F.isOne( alpha) &&
					    !F.isMOne( alpha))
						// Fix-up: compute C *= alpha
						fscalin(F,m,n,alpha,C,ldc);
				}
			       	else {
					WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha,
						  A, lda, B, ldb, beta, C, ldc, kmax,w,base);
					DynamicPealing (F, ta, tb, m, n, k, alpha,
							A, lda, B, ldb, beta, C, ldc, kmax);
				}
			}
		}


		template <>
		inline void WinoMain (const FFPACK:: Modular<float>& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t m, const size_t n, const size_t k,
				      const float alpha,
				      const float* A, const size_t lda,
				      const float* B, const size_t ldb,
				      const float beta,
				      float * C, const size_t ldc,
				      const size_t kmax, const size_t w, const FFLAS_BASE base)
		{
			if (w == 0) {
				ClassicMatmul (F, ta, tb, m, n, k, alpha,
					       A, lda, B, ldb, beta, C, ldc, kmax,base);
			}
			else {
				if (k <= kmax) { // switch on float
					FloatDomain::Element _alpha, _beta;
					_beta = beta;
					if (F.isMOne( alpha)) _alpha = -1.0;
					else {
						// Compute C = A*B + beta/alpha.C
						// and then C *= alpha
						if (! F.isOne( alpha)) {
							FFLASFFPACK_check(!F.isZero(alpha));
							F.divin (_beta, alpha);
						}
						_alpha = 1.0;
					}
#if 0 /* BB: inutile */
					size_t ma, ka, kb, nb; //mb, na;
					if (ta == FflasTrans) { ma = k; ka = m; }
					else { ma = m; ka = k; }
					if (tb == FflasTrans) { kb = n; nb = k; }
					else {  kb = k; nb = n; }
#endif
					// recursive call
					WinoMain (FloatDomain(), ta, tb, m, n, k, _alpha,
						  A, lda, B, ldb, _beta, C, ldc, kmax, w,base);
					// Conversion float = >  GFq
					finit(F,m,n,C,ldc);

					if (!F.isOne( alpha) &&
					    !F.isMOne( alpha))
						// Fix-up: compute C *= alpha
						fscalin(F,m,n,alpha,C,ldc);
				}
			       	else{
					WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha,
						  A, lda, B, ldb, beta, C, ldc, kmax,w,base);
					DynamicPealing (F, ta, tb, m, n, k, alpha,
							A, lda, B, ldb, beta, C, ldc, kmax);
				}
			}
		}

		template  < class Field >
		inline void
		DynamicPealing (const Field& F,
				const FFLAS_TRANSPOSE ta,
				const FFLAS_TRANSPOSE tb,
				const size_t m, const size_t n, const size_t k,
				const typename Field::Element alpha,
				const typename Field::Element* A, const size_t lda,
				const typename Field::Element* B, const size_t ldb,
				const typename Field::Element beta,
				typename Field::Element* C, const size_t ldc,
				const size_t kmax)
		{
			const typename Field::Element *a12, *a21, *b12, *b21;
			size_t inca12, inca21, incb12, incb21, ma, na, mb, nb;
			size_t mkn = (n & 0x1)+ ((k & 0x1) << 1)+  ((m & 0x1) << 2);

			if (ta == FflasTrans) {
				ma = k;
				na = m;
				a12 = A+(k-1)*lda;
				inca12 = 1;
				a21 = A+m-1;
				inca21 = lda;
			}
			else {
				ma = m;
				na = k;
				a12 = A+k-1;
				inca12 = lda;
				a21 = A+(m-1)*lda;
				inca21 = 1;
			}
			if (tb == FflasTrans) {
				mb = n;
				nb = k;
				b12 = B+(n-1)*ldb;
				incb12 = 1;
				b21 = B+k-1;
				incb21 = ldb;
			}
			else {
				mb = k;
				nb = n;
				b12 = B+n-1;
				incb12 = ldb;
				b21 = B+(k-1)*ldb;
				incb21 = 1;
			}
			switch (mkn) {
			case 1: // n oddsized
				fgemv (F, ta, ma, na, alpha, A, lda, b12, incb12, beta, C+n-1,ldc);
				break;

			case 2: // k oddsized
				fger (F, m, n, alpha, a12, inca12, b21, incb21, C, ldc);
				break;

			case 3: // n, k oddsized
				fgemv (F, ta, ma, na, alpha, A, lda, b12, incb12, beta, C+n-1,ldc);
				fger (F, m, n-1, alpha, a12, inca12, b21, incb21, C, ldc);
				break;

			case 4: // m oddsized
				fgemv(F, (tb == FflasTrans)?FflasNoTrans:FflasTrans, mb, nb,
				      alpha, B, ldb, a21, inca21, beta, C+(m-1)*ldc, 1);
				break;

			case 5: // m, n oddsized
				if (tb == FflasTrans)
					mb--;
				else
					nb--;
				fgemv (F, ta, ma, na, alpha, A, lda, b12, incb12, beta, C+n-1, ldc);
				fgemv (F, (tb==FflasTrans)?FflasNoTrans:FflasTrans, mb, nb,
				       alpha, B, ldb, a21, inca21, beta, C+(m-1)*ldc, 1);
				break;

			case 6: // m, k oddsized
				fger (F, m-1, n, alpha, a12, inca12, b21, incb21, C, ldc);
				fgemv(F, (tb==FflasTrans)?FflasNoTrans:FflasTrans, mb, nb,
				      alpha, B, ldb, a21, inca21, beta, C+(m-1)*ldc, 1);
				break;

			case 7: // m, k, n oddsized
				if (tb == FflasTrans)
					mb--;
				else
					nb--;
				// Block NW
				fger (F, m-1, n-1, alpha, a12, inca12, b21, incb21, C, ldc);
				// Block SW
				fgemv (F, (tb==FflasTrans)?FflasNoTrans:FflasTrans, mb, nb,
				       alpha, B, ldb, a21, inca21, beta, C+(m-1)*ldc, 1);
				// Block NE
				fgemv (F, ta, ma, na, alpha, A, lda, b12, incb12, beta, C+n-1, ldc);
				break;
			}
		}

	} // Protected Winomain

	// Unsafe matmul over Z
	// For internal usage only (or use it with care)
	template<>
	inline double*
	fgemm< FFPACK:: UnparametricField<double> > ( const  FFPACK:: UnparametricField<double>& F,
					    const FFLAS_TRANSPOSE ta,
					    const FFLAS_TRANSPOSE tb,
					    const size_t m,
					    const size_t n,
					    const size_t k,
					    const double alpha,
					    const double* A, const size_t lda,
					    const double* B, const size_t ldb,
					    const double beta,
					    double* C, const size_t ldc,
					    const size_t w)
	{

		if (!(m && n && k)) return C;

		if (F.isZero (alpha)){
			fscalin(F,m,n,beta,C,ldc);
			return C;
		}


		Protected::WinoMain (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta,
			  C, ldc, k+1, w, FflasDouble);
		return C;
	}

	template<>
	inline float*
	fgemm< FFPACK:: UnparametricField<float> > ( const  FFPACK:: UnparametricField<float>& F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m,
					   const size_t n,
					   const size_t k,
					   const float alpha,
					   const float* A, const size_t lda,
					   const float* B, const size_t ldb,
					   const float beta,
					   float* C, const size_t ldc,
					   const size_t w)
	{

		if (!(m && n && k)) return C;

		if (F.isZero (alpha)){
			fscalin(F,m,n,beta,C,ldc);
			return C;
		}

		Protected::WinoMain (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta,
			  C, ldc, k+1, w, FflasFloat);
		return C;
	}

	template<>
	inline double*
	fgemm< FFPACK:: UnparametricField<double> > (const  FFPACK:: UnparametricField<double>& F,
					   const FFLAS_TRANSPOSE ta,
					   const FFLAS_TRANSPOSE tb,
					   const size_t m,
					   const size_t n,
					   const size_t k,
					   const double alpha,
					   const double* A, const size_t lda,
					   const double* B, const size_t ldb,
					   const double beta,
					   double* C, const size_t ldc)
	{
		return fgemm (F, ta, tb, m, n ,k, alpha, A, lda, B, ldb, beta, C, ldc, WinoSteps (std::min(m,std::min(k,n))));
	}

	template<>
	inline float*
	fgemm< FFPACK:: UnparametricField<float> > (const  FFPACK:: UnparametricField<float>& F,
					  const FFLAS_TRANSPOSE ta,
					  const FFLAS_TRANSPOSE tb,
					  const size_t m,
					  const size_t n,
					  const size_t k,
					  const float alpha,
					  const float* A, const size_t lda,
					  const float* B, const size_t ldb,
					  const float beta,
					  float* C, const size_t ldc)
	{
		return fgemm (F, ta, tb, m, n ,k, alpha, A, lda, B, ldb, beta, C, ldc, WinoSteps (std::min(m,std::min(k,n))));
	}


	template < class Field >
	inline typename Field::Element*
	fsquare (const Field& F,
		 const FFLAS_TRANSPOSE ta,
		 const size_t n, const typename Field::Element alpha,
		 const typename Field::Element* A, const size_t lda,
		 const typename Field::Element beta,
		 typename Field::Element* C, const size_t ldc)
	{

		double alphad, betad;
		F.convert (alphad, alpha);
		if (F.isMOne (beta))
			betad = -1.0;
		else
			F.convert (betad, beta);

		// Double  matrices initialisation
		DoubleDomain::Element * Ad = new DoubleDomain::Element[n*n];
		DoubleDomain::Element * Cd = new DoubleDomain::Element[n*n];
		// Conversion finite Field = >  double
		Protected::MatF2MatD (F, Ad, n, A, lda, n, n);
		if (!F.isZero(beta)) Protected::MatF2MatD (F, Cd, n, C, ldc, n, n);

		// Call to the blas Multiplication
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
			     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
			     (DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
			     (DoubleDomain::Element) betad, Cd, (int)n);
		// Conversion double = >  Finite Field
		delete[] Ad;
		Protected::MatD2MatF (F, C, ldc, Cd, n, n, n);
		delete[] Cd;
		return C;
	}

	template <>
	inline double* fsquare (const  FFPACK:: ModularBalanced<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		if (C==A) {
			double * Ad = new double[n*n];
			for (size_t i=0; i < n; ++i)
				fcopy (F, n,Ad+i*n, 1, A+i*lda, 1);
			fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
			delete[] Ad;
		}
	       	else
			fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
		// Conversion double = >  Finite Field
		finit(F,n,n,C,ldc);
		return C;
	}

	template <>
	inline float * fsquare (const  FFPACK:: ModularBalanced<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		if (C==A) {
			float * Ad = new float[n*n];
			for (size_t i=0; i < n; ++i)
				fcopy (F, n,Ad+i*n, 1, A+i*lda, 1);
			fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
			delete[] Ad;
		}
	       	else
			fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
		// Conversion float = >  Finite Field
		finit(F,n,n,C,ldc);
		return C;
	}

	template <>
	inline double* fsquare (const  FFPACK:: Modular<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		if (C==A) {
			double * Ad = new double[n*n];
			for (size_t i=0; i < n; ++i)
				fcopy (F, n,Ad+i*n, 1, A+i*lda, 1);
			fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
			delete[] Ad;
		}
	       	else
			fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
		// Conversion double = >  Finite Field
		finit(F,n,n,C,ldc);
		return C;
	}

	template <>
	inline float * fsquare (const  FFPACK:: Modular<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		if (C==A) {
			float * Ad = new float[n*n];
			for (size_t i=0; i < n; ++i)
				fcopy (F, n,Ad+i*n, 1, A+i*lda, 1);
			fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
			delete[] Ad;
		}
	       	else
			fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
		// Conversion float = >  Finite Field
		finit(F,n,n,C,ldc);
		return C;
	}

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
