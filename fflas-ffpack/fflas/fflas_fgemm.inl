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
			else{
#ifdef OLD_DYNAMIC_PEALING
				WinoCalc (F, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb,
					  beta, C, ldc, kmax, w,base);
				DynamicPealing (F, ta, tb, m, n, k, alpha, A, lda, B, ldb,
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
				// std::cout  << "w =" << ww << "m =" << m << ", mr = " << mr  << ", m2 = " << m2 << std::endl;

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
			size_t mkn = (bool)(nr > 0)+ ((bool)(kr > 0) << 1)+  ((bool)(mr > 0) << 2);
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
