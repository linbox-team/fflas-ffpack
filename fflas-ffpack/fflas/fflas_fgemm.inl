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

// #include "fflas_fgemm/matmul_algos.inl"
#include "fflas_fgemm/winograd.inl"
#include "fflas_fgemm/winograd_acc.inl"
#include "fflas_fgemm/winograd_acc2.inl"
#include "fflas_fgemm/winograd_ip.inl"
// #include "fflas_fgemm/bini.inl"

#define NEWWINO
// #define OLDWINO

namespace FFLAS {

	DoubleDomain associatedDomain (const FFPACK::Modular<double> & )
	{
		return DoubleDomain();
	}

	DoubleDomain associatedDomain (const FFPACK::ModularBalanced<double> & )
	{
		return DoubleDomain();
	}

	FloatDomain associatedDomain (const FFPACK::Modular<float> & )
	{
		return FloatDomain();
	}

	FloatDomain associatedDomain (const FFPACK::ModularBalanced<float> & )
	{
		return FloatDomain();
	}

	namespace Protected {

		// G is (Float/Double)Domain
		template  < class Field, class FloatField >
		inline void ClassicMatmulFloat (const Field& F,
						const FFLAS_TRANSPOSE ta,
						const FFLAS_TRANSPOSE tb,
						const size_t m, const size_t n,const size_t k,
						const typename Field::Element alpha,
						const typename Field::Element * A, const size_t lda,
						const typename Field::Element * B, const size_t ldb,
						const typename Field::Element beta,
						typename Field::Element* C, const size_t ldc,
						const size_t kmax, const FFLAS_BASE base
						, const FloatField & G, const size_t k2)
		{
			typedef typename FloatField::Element FloatElement ;
			FloatElement alphad, betad;
			FloatElement * Add = new FloatElement[m*k2];
			FloatElement * Bdd = new FloatElement[k2*n];
			FloatElement * Cd = new FloatElement[m*n];

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
					typename Field::Element tmp;
					FFLASFFPACK_check(!F.isZero(alpha));
					F.div (tmp, beta, alpha);
					F.convert (betad, tmp);
				}
			}

			size_t dlda, dldb;
			if (!F.isZero(beta))
				fconvert (F, m, n, Cd, n, C, ldc);

			if (ta == FflasTrans) {
				dlda = m;
				fconvert(F, remblock, m, Add, dlda, A+k2*nblock*lda, lda);
			}
			else {
				dlda = k2;
				fconvert(F, m, remblock, Add, dlda, A+k2*nblock, lda);
			}
			if (tb == FflasTrans) {
				dldb = k2;
				fconvert (F, n, remblock, Bdd, k2, B+k2*nblock, ldb);
			}
			else {
				dldb = n;
				fconvert(F, remblock, n, Bdd, dldb, B+k2*nblock*ldb, ldb);
			}

			ClassicMatmul (G, ta, tb, m, n, remblock, alphad, Add, dlda,
				       Bdd, dldb, betad, Cd, n, kmax,base );

			finit (F, m, n, C, ldc, Cd, n);
			fconvert(F, m, n, Cd, n, C, ldc);

			for (size_t i = 0; i < nblock; ++i) {
				if (ta == FflasTrans) fconvert(F, k2, m, Add, dlda, A+k2*i*lda, lda);
				else fconvert(F, m, k2, Add, dlda,  A+k2*i, lda);
				if (tb == FflasTrans) fconvert(F, n, k2, Bdd, dldb, B+k2*i, ldb);
				else fconvert(F, k2, n, Bdd, dldb, B+k2*i*ldb, ldb);

				ClassicMatmul (G, ta, tb, m, n, k2, alphad, Add, dlda,
					       Bdd, dldb, 1.0, Cd, n, kmax,base);
				finit(F, m, n, C, ldc, Cd, n);
				fconvert(F, m, n, Cd, n, C, ldc);
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
			delete[] Add;
			delete[] Bdd;
			delete[] Cd;

		}


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

			size_t k2 = std::min(k,kmax); // Size of the blocks

			if (k2 > 1) {
				if (base == FflasFloat) {
					FloatDomain G ;
					ClassicMatmulFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base,G,k2);
				}
				else { // base == FflasDouble
					DoubleDomain G ;
					ClassicMatmulFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base,G,k2);
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

			FFLASFFPACK_check(lda);
			FFLASFFPACK_check(ldb);
			FFLASFFPACK_check(ldc);
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
			FFLASFFPACK_check(lda);
			FFLASFFPACK_check(ldb);
			FFLASFFPACK_check(ldc);

			cblas_sgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
				     (int)m, (int)n, (int)k, (FloatDomain::Element) alpha,
				     Ad, (int)lda, Bd, (int)ldb, (FloatDomain::Element) beta,Cd, (int)ldc);
		}


		// F is Modular(Balanced)<float/double>
		template<class Field>
		inline void ClassicMatmulCommon (const Field & F,
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
			typename Field::Element _alpha, _beta;
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

			ClassicMatmul (associatedDomain(F), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
				       B+nblock*shiftB, ldb, _beta, C, ldc, kmax,base);
			finit(F,m,n,C,ldc);
			for (size_t i = 0; i < nblock; ++i) {
				ClassicMatmul (associatedDomain(F), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
					       B+i*shiftB, ldb, F.one, C, ldc, kmax,base);
				finit(F,m,n,C,ldc);
			}
			if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
				fscalin(F,m,n,alpha,C,ldc);
			}
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
			return ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base);
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
			return ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base);
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
			return ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base);
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
			return ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,base);
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

			if (F.isZero(beta)) {
				if (kr == nr && ta == FflasNoTrans && tb == FflasNoTrans) {

					std::cout << "tested" << std::endl;
					typedef typename Field::Element Element ;
					size_t ldA = lda;
					size_t ldB = ldb;
					Element * Ac = new Element[mr*2*ldA] ;
					Element * Bc = new Element[kr*2*ldB] ;
					fcopy(F,mr*2,kr*2,Ac,ldA,A,lda);
					fcopy(F,kr*2,nr*2,Bc,ldB,B,ldb);

					BLAS3::Winograd(F,ta,tb,mr,nr,kr,alpha,Ac,ldA,Bc,ldB,beta,C,ldc,kmax,w,base);
					delete[] Ac;
					delete[] Bc;
					std::cout << "done" << std::endl;
				}
				else
				{
					BLAS3::Winograd(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
				}
			}
			else
#ifdef NEWWINO
				BLAS3::WinogradAcc(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
#elif defined(OLDWINO)
				BLAS3::WinogradAccOld(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
#elif defined(NEWACC)
				// BLAS3::WinogradAcc2(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
				BLAS3::WinogradAcc3(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
#else
#error "you need to make a choice for a BLAS3 mat mul schedule"
#endif

		}

#define OLD_DYNAMIC_PEALING
		// dispatches according to w = 0 or not
		template<class Field>
		inline  void WinoMainGeneric (const Field& F,
				       const FFLAS_TRANSPOSE ta,
				       const FFLAS_TRANSPOSE tb,
				       const size_t m, const size_t n, const size_t k,
				       const typename Field::Element alpha,
				       const typename Field::Element * A, const size_t lda,
				       const typename Field::Element * B, const size_t ldb,
				       const typename Field::Element beta,
				       typename Field::Element * C, const size_t ldc,
				       const size_t kmax, const size_t w, const FFLAS_BASE base)
		{

			if (!m || !n )
				return;
			if (!k)
				return fscalin(F,m,n,beta,C,ldc);

			if (w == 0) {
				ClassicMatmul (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
					       beta, C, ldc, kmax,base);
			}
			else {
#ifdef OLD_DYNAMIC_PEALING
				WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb,
					  beta, C, ldc, kmax, w,base);

				FFLASFFPACK_check(m-(m/2)*2 == (m&0x1));
				FFLASFFPACK_check(n-(n/2)*2 == (n&0x1));
				FFLASFFPACK_check(k-(k/2)*2 == (k&0x1));

				DynamicPealing (F, ta, tb, m, n, k, m&0x1, n&0x1, k&0x1, alpha, A, lda, B, ldb,
						beta, C, ldc, kmax);
#else
				size_t ww = w ;
				size_t m2 = (m >> ww) << (ww-1) ;
				size_t n2 = (n >> ww) << (ww-1) ;
				size_t k2 = (k >> ww) << (ww-1) ;

				WinoCalc (F, ta, tb, m2, n2, k2, alpha, A, lda, B, ldb,
					  beta, C, ldc, kmax, w,base);

				size_t mr = m -2*m2;
				size_t nr = n -2*n2;
				size_t kr = k -2*k2;

				FFLASFFPACK_check(m == m2*2+mr);
				FFLASFFPACK_check(n == n2*2+nr);
				FFLASFFPACK_check(k == k2*2+kr);

				DynamicPealing2 (F, ta, tb, m, n, k, mr, nr, kr, alpha, A, lda, B, ldb,
						beta, C, ldc, kmax);


#endif
			}
		}

		// G is (Float/Double)Domain
		template <class Field, class FloatField>
		inline void  WinoMainFloat (const Field& F,
					    const FFLAS_TRANSPOSE ta,
					    const FFLAS_TRANSPOSE tb,
					    const size_t m, const size_t n, const size_t k,
					    const typename Field::Element alpha,
					    const typename Field::Element* A,const size_t lda,
					    const typename Field::Element* B,const size_t ldb,
					    const typename Field::Element beta,
					    typename Field::Element * C, const size_t ldc,
					    const size_t kmax, const size_t w,
					    const FFLAS_BASE base
					    , const FloatField &G)
		{
			FFLASFFPACK_check(lda);
			FFLASFFPACK_check(ldb);
			FFLASFFPACK_check(ldc);
			typedef typename FloatField::Element FloatElement;
			FloatElement alphad, betad;
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
			FloatElement * Ad = new FloatElement[m*k];
			FloatElement * Bd = new FloatElement[k*n];
			FloatElement * Cd = new FloatElement[m*n];
			// Conversion GFq = >  double
			size_t ma, ka, kb, nb; //mb, na
			if (ta == FflasTrans) { ma = k; ka = m; }
			else { ma = m; ka = k; }
			if (tb == FflasTrans) { kb = n; nb = k; }
			else {  kb = k; nb = n; }

			fconvert(F, ma, ka, Ad, ka, A, lda);
			fconvert(F, kb, nb, Bd, nb, B, ldb);
			if (!F.isZero(beta))
				fconvert(F, m, n, Cd, n, C, ldc);
			// recursive call
			WinoMain (G, ta, tb, m, n, k, alphad,
				  Ad, ka, Bd, nb, betad, Cd, n, kmax, w,base);
			// Conversion double = >  GFq
			finit(F, m, n, C, ldc, Cd, n);

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

		// Control the switch with classic multiplication
		// Fix-up for odd-sized matrices using dynamic pealing
		// for matrices over double
		template <>
		inline  void WinoMain (const DoubleDomain& F,
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
			WinoMainGeneric(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
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
			WinoMainGeneric(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
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
			if (w > 0 && k <= kmax) {
				if (base == FflasDouble){
					DoubleDomain G ;
					return WinoMainFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base,G);
				}
				else { // FloatDomain
					FloatDomain G ;
					return WinoMainFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base,G);
				}
			}
			// w = 0 or k> kmax
			WinoMainGeneric(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);

		}


		// F is Modular(Balanced)<float/double>
		template <class Field>
		inline void WinoMainCommon (const Field& F,
					    const FFLAS_TRANSPOSE ta,
					    const FFLAS_TRANSPOSE tb,
					    const size_t m, const size_t n, const size_t k,
					    const typename Field::Element alpha,
					    const typename Field::Element* A, const size_t lda,
					    const typename Field::Element* B, const size_t ldb,
					    const typename Field::Element beta,
					    typename Field::Element * C, const size_t ldc,
					    const size_t kmax, const size_t w, const FFLAS_BASE base)
		{
			if (w > 0 && k <= kmax) { // switch on delayed modulus
				typename Field::Element _alpha, _beta;
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
				// recursive call
				WinoMain (associatedDomain(F), ta, tb, m, n, k, _alpha,
					  A, lda, B, ldb, _beta, C, ldc, kmax, w,base);
				// Modular reduction
				finit(F,m,n,C,ldc);
				if (!F.isOne( alpha ) && !F.isMOne( alpha ))
					// Fix-up: compute C *= alpha
					fscalin(F,m,n,alpha,C,ldc);

				return;
			}
			// w = 0 or k> kmax
			WinoMainGeneric(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);

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
			WinoMainCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
		}

		template <>
		inline void WinoMain (const FFPACK:: ModularBalanced<float>& F,
				      const FFLAS_TRANSPOSE ta,
				      const FFLAS_TRANSPOSE tb,
				      const size_t m, const size_t n, const size_t k,
				      const float alpha,
				      const float * A,const size_t lda,
				      const float * B,const size_t ldb,
				      const float beta,
				      float * C, const size_t ldc,
				      const size_t kmax, const size_t w,
				      const FFLAS_BASE base)
		{
			WinoMainCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
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
			WinoMainCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
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
			WinoMainCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,kmax,w,base);
		}


		template  < class Field >
		inline void
		DynamicPealing (const Field& F,
				const FFLAS_TRANSPOSE ta,
				const FFLAS_TRANSPOSE tb,
				const size_t m, const size_t n, const size_t k,
				const size_t mr, const size_t nr, const size_t kr,
				const typename Field::Element alpha,
				const typename Field::Element* A, const size_t lda,
				const typename Field::Element* B, const size_t ldb,
				const typename Field::Element beta,
				typename Field::Element* C, const size_t ldc,
				const size_t kmax)
		{
			const typename Field::Element *a12, *a21, *b12, *b21;
			size_t inca12, inca21, incb12, incb21, ma, na, mb, nb;
			size_t mkn = nr + (kr << 1)+  (mr << 2);

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

		template  < class Field >
		inline void
		DynamicPealing2 (const Field& F,
				 const FFLAS_TRANSPOSE ta,
				 const FFLAS_TRANSPOSE tb,
				 const size_t m, const size_t n, const size_t k,
				 const size_t mr, const size_t nr, const size_t kr,
				 const typename Field::Element alpha,
				 const typename Field::Element* A, const size_t lda,
				 const typename Field::Element* B, const size_t ldb,
				 const typename Field::Element beta,
				 typename Field::Element* C, const size_t ldc,
				 const size_t kmax)
		{
			size_t mkn =(size_t)( (bool)(nr > 0)+ ((bool)(kr > 0) << 1)+  ((bool)(mr > 0) << 2));
			if (mkn == 0) return;

			const typename Field::Element *a12, *a21, *b12, *b21;
			if (ta == FflasTrans) {
				a12 = A+(k-kr)*lda;
				a21 = A+(m-mr);
			}
			else {
				a12 = A+(k-kr);
				a21 = A+(m-mr)*lda;
			}
			if (tb == FflasTrans) {
				b12 = B+(n-nr)*ldb;
				b21 = B+(k-kr);
			}
			else {
				b12 = B+(n-nr);
				b21 = B+(k-kr)*ldb;
			}
			switch (mkn) {
			case 1: // n oddsized
				fgemm (F, ta, tb, m, nr, k, alpha, A, lda, b12, ldb, beta, C+(n-nr), ldc);

				break;

			case 2: // k oddsized
				fgemm (F, ta, tb, m, n, kr, alpha, a12, lda, b21, ldb, F.one, C, ldc);

				break;

			case 3: // n, k oddsized
				fgemm (F, ta, tb, m, nr, k, alpha, A, lda, b12, ldb, beta, C+(n-nr), ldc);
				fgemm (F, ta, tb, m, n-nr, kr, alpha, a12, lda, b21, ldb, F.one, C, ldc);

				break;

			case 4: // m oddsized
				fgemm (F,  ta, tb, mr, n, k, alpha, a21, lda, B, ldb, beta, C+(m-mr)*ldc, ldc);

				break;

			case 5: // m, n oddsized
				fgemm (F, ta, tb, m, nr, k, alpha, A, lda, b12, ldb, beta, C+(n-nr), ldc);
				fgemm (F, ta, tb, mr, n-nr, k, alpha, a21, lda, B, ldb, beta, C+(m-mr)*ldc, ldc);

				break;

			case 6: // m, k oddsized
				fgemm (F, ta, tb, m-mr, n, kr, alpha, a12, lda, b21, ldb, F.one, C, ldc);
				fgemm (F,  ta, tb, mr, n, k, alpha, a21, lda, B, ldb, beta, C+(m-mr)*ldc, ldc);

				break;

			case 7: // m, k, n oddsized
				// Block NW
				fgemm (F, ta, tb, m-mr, n-nr, kr, alpha, a12, lda, b21, ldb, F.one, C, ldc);
				// Block SW
				fgemm (F,  ta, tb, mr, n-nr, k, alpha, a21, lda, B, ldb, beta, C+(m-mr)*ldc, ldc);
				// Block NE
				fgemm (F, ta, tb, m, nr, k, alpha, A, lda, b12, ldb, beta, C+(n-nr), ldc);

				break;
			}
		}


	} // Protected Winomain

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
		fconvert (F, n, n, Ad, n, A, lda);
		if (!F.isZero(beta)) fconvert(F, n, n, Cd, n, C, ldc);

		// Call to the blas Multiplication
		FFLASFFPACK_check(n);
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
			     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
			     (DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
			     (DoubleDomain::Element) betad, Cd, (int)n);
		// Conversion double = >  Finite Field
		delete[] Ad;
		finit (F,n,n, C, ldc, Cd, n);
		delete[] Cd;
		return C;
	}

	namespace Protected {

		// F is Modular(Balanced)<float/double>
		template < class Field >
		inline typename Field::Element*
		fsquareCommon (const Field& F,
			       const FFLAS_TRANSPOSE ta,
			       const size_t n, const typename Field::Element alpha,
			       const typename Field::Element* A, const size_t lda,
			       const typename Field::Element beta,
			       typename Field::Element* C, const size_t ldc)
		{
			typedef typename Field::Element Element ; // double or float
			if (C==A) {
				Element * Ad = new Element[n*n];
				fcopy(F,n,n,Ad,n,A,lda);
				fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
				delete[] Ad;
			}
			else
				fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
			// Conversion double/float = >  Finite Field
			finit(F,n,n,C,ldc);
			return C;

		}

	} // Protected

	template <>
	inline double* fsquare (const  FFPACK:: ModularBalanced<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: ModularBalanced<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline double* fsquare (const  FFPACK:: Modular<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: Modular<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
