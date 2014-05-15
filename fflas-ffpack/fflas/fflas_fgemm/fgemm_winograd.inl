/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
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

/** @file fflas_fgemm/fgemm_winograd.h
 * @brief Strassen--Winograd matrix multiplication.
 * @warning The domain is supposed to be a field since some divisions are required for efficiency purposes
 * An alternative has to be written for finite rings if necessary
 */

#ifndef __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL
#define __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL

#include "fflas_bounds_winograd.inl"
#include "fgemm_classical.inl"

namespace FFLAS {

	template<typename FieldTrait>
	struct Winograd2Helper : public MMParameters {
		int w ;
		size_t kmax ;
		FFLAS_BASE base ;
		// double pmin,pmax ;
		// bool   delay ;

		Winograd2Helper() :
			w(-1) , kmax(0), base(FflasDouble)
			// ,delay(true),pmin(0),pmax(0)
		{}

		Winograd2Helper(int rec) :
			w(rec) , kmax(0), base(FflasDouble)
			// ,delay(true),pmin(0),pmax(0)
		{}
		template<class FT2>
		Winograd2Helper(Winograd2Helper<FT2> WH) :
			w(WH.w), kmax(WH.kmax), base(WH.base)
		{}

		// copy constructor
		// compute parameters (full delay)
		template<class Field>
		void computeParameters(const Field & F
				       , const size_t & m
				       , const size_t & n
				       , const size_t & k
				       , const typename Field::Element &alpha
				       , const typename Field::Element &beta)
		{
			// if (delay) {
			bool winoLevelProvided = (w != (int(-1)));
			// size_t kmax = 0;
			int winolevel = w;
			// FFLAS_BASE base;
			typename Field::Element gamma;
			F.div(gamma,beta,alpha);
			Protected::MatMulParametersWinograd (F, m, n, k, gamma, kmax, base,
							     winolevel, winoLevelProvided);
			w = winolevel;
			// }
			// else {
			// compute kmax
			// pmin = min(F);
			// pmax = max(F);
			// }

		}

#if 0
		// tell if one can do some winograd step
		template<class Field>
		bool delay(const Field & F
			   , const size_t & m
			   , const size_t & n
			   , const size_t & k
			   , const typename Field::Element &alpha
			   , const typename Field::Element &beta)
		;
#endif





	};
} // FFLAS

#include "schedule_winograd.inl"
#include "schedule_winograd_acc.inl"
#include "schedule_winograd_acc_ip.inl"
#include "schedule_winograd_ip.inl"
// #include "fflas_fgemm/bini.inl"

#ifndef NEWWINO
#define NEWWINO
#endif

// #define OLDWINO


// DynamicPealing,  WinogradMain, WinogradCalc
namespace FFLAS { namespace Protected {

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

	// #define NEWIP
	// #define NEWACCIP

	// Switch between the scheduling for Strassen-Winograd Multiplication
	template < class Field >
	inline void WinogradCalc (const Field& F,
				  const FFLAS_TRANSPOSE ta,
				  const FFLAS_TRANSPOSE tb,
				  const size_t mr, const size_t nr, const size_t kr,
				  const typename Field::Element alpha,
				  const typename Field::Element* A,const size_t lda,
				  const typename Field::Element* B,const size_t ldb,
				  const typename Field::Element beta,
				  typename Field::Element * C, const size_t ldc,
				  const Winograd2Helper<typename FieldTraits<Field>::value> & H)
	{

#if defined(NEWIP) or defined(NEWACCIP)  /*  XXX TESTS ONLY */
		typedef typename Field::Element Element ;
		Element * Ac;
		Element * Bc;
		if (ta == FflasNoTrans) {
			Ac = new Element[mr*2*lda] ;
			fcopy(F,mr*2,kr*2,A,lda,Ac,lda);
		}
		else {
			Ac = new Element[kr*2*lda] ;
			fcopy(F,kr*2,mr*2,A,lda,Ac,lda);
		}
		if (tb == FflasNoTrans) {
			Bc = new Element[kr*2*ldb] ;
			fcopy(F,kr*2,nr*2,B,ldb,Bc,ldb);
		}
		else {
			Bc = new Element[nr*2*ldb] ;
			fcopy(F,nr*2,kr*2,B,ldb,Bc,ldb);
		}
#endif

		if (F.isZero(beta)) {
#ifdef NEWIP /*  NOT IP --- TESTS ONLY */
			// (kr == nr  && kr <= mr /*  if not transposed */)
			// we copy because they erase stuff
			// bool normal =  (ta == FflasNoTrans && tb == FflasNoTrans) ;
			bool normal = true;

			// std::cout << (ta==FflasNoTrans) << ',' << (tb==FflasNoTrans) << std::endl;

			if (kr == nr && kr == mr && normal) {
				// BLAS3::Winograd_L_S(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
				// BLAS3::Winograd_R_S(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
				BLAS3::Winograd_LR_S(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
			}
			else
#endif
			{
				BLAS3::Winograd(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,H);
			}

		}
		else {
#ifdef NEWACCIP /*  test only */
			// std::cout << (ta==FflasNoTrans) << ',' << (tb==FflasNoTrans) << std::endl;
			if (kr == nr && kr == mr ) {
				// std::cout << 'h' << std::endl;
				BLAS3::WinogradAcc_L_S(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
				// BLAS3::WinogradAcc_R_S(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
			}
			else {
				BLAS3::WinogradAcc_LR(F,ta,tb,mr,nr,kr,alpha,Ac,lda,Bc,ldb,beta,C,ldc,H);
			}


#elif defined(NEWWINO)
			BLAS3::WinogradAcc_3_21(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,H);
#elif defined(OLDWINO)
			BLAS3::WinogradAcc_3_23(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,H);
#elif defined(NEWACC)
			// BLAS3::WinogradAcc_2_24(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,H);
			BLAS3::WinogradAcc_2_27(F,ta,tb,mr,nr,kr,alpha,A,lda,B,ldb,beta,C,ldc,H);
#else
#error "you need to make a choice for a BLAS3 mat mul schedule"
#endif

		}
#if defined(NEWIP) or defined(NEWACCIP)  /*  NOT IP --- TESTS ONLY */
		delete[] Ac;
		delete[] Bc;
#endif

	} // WinogradCalc

	//#define OLD_DYNAMIC_PEALING

}// namespace Protected
} // FFLAS

namespace FFLAS {

template<class Field>
inline  void fgemm2 (const Field& F,
		     const FFLAS_TRANSPOSE ta,
		     const FFLAS_TRANSPOSE tb,
		     const size_t m, const size_t n, const size_t k,
		     const typename Field::Element alpha,
		     const typename Field::Element * A, const size_t lda,
		     const typename Field::Element * B, const size_t ldb,
		     const typename Field::Element beta,
		     typename Field::Element * C, const size_t ldc,
		     const Winograd2Helper<typename FieldTraits<Field>::value> & H)
{
	if (!m || !n ) return;

	if (!k)	return fscalin(F,m,n,beta,C,ldc);

	if (H.w == 0) return fgemm2(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, ClassicHelper<typename FieldTraits<Field>::value>(H.kmax,H.base));

	// Then w >0
	if (Protected::AreEqual< typename FieldTraits<Field>::value,
	    FieldCategories::FloatingPointConvertibleTag>::value){
		// Field is convertible to a floating point representation
		if (k <= H.kmax) {
			if (H.base == FflasDouble)
				return Protected::fgemm_convert<double,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
			else // FloatDomain
				return Protected::fgemm_convert<float,Field>(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
		}

	} else if (Protected::AreEqual<typename FieldTraits<Field>::value,
		   FieldCategories::ModularFloatingPointTag>::value){
		// Field is a Modular[Balanced]<float,double>

		if (k <= H.kmax) { // switch on delayed modulus
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
			// call on Z with no modulo
			fgemm2 (associatedDomain(F), ta, tb, m, n, k, _alpha,
				A, lda, B, ldb, _beta, C, ldc,
				Winograd2Helper<FieldCategories::FloatingPointTag>(H));
			// Modular reduction
			finit(F,m,n,C,ldc);
			if (!F.isOne( alpha ) && !F.isMOne( alpha ))
				// Fix-up: compute C *= alpha
				fscalin(F,m,n,alpha,C,ldc);
			return;
		}
	}

#ifdef OLD_DYNAMIC_PEALING
	WinogradCalc (F, ta, tb, m/2, n/2, k/2, alpha, A, lda, B, ldb, beta, C, ldc, H);

	FFLASFFPACK_check(m-(m/2)*2 == (m&0x1));
	FFLASFFPACK_check(n-(n/2)*2 == (n&0x1));
	FFLASFFPACK_check(k-(k/2)*2 == (k&0x1));

	Protected::DynamicPealing (F, ta, tb, m, n, k, m&0x1, n&0x1, k&0x1, alpha, A, lda, B, ldb,
				   beta, C, ldc, H.kmax);
#else
	size_t ww = H.w ;
	size_t m2 = (m >> ww) << (ww-1) ;
	size_t n2 = (n >> ww) << (ww-1) ;
	size_t k2 = (k >> ww) << (ww-1) ;

	Protected::WinogradCalc (F, ta, tb, m2, n2, k2, alpha, A, lda, B, ldb,
				 beta, C, ldc, H);

	size_t mr = m -2*m2;
	size_t nr = n -2*n2;
	size_t kr = k -2*k2;

	FFLASFFPACK_check(m == m2*2+mr);
	FFLASFFPACK_check(n == n2*2+nr);
	FFLASFFPACK_check(k == k2*2+kr);

	Protected::DynamicPealing2 (F, ta, tb, m, n, k, mr, nr, kr, alpha, A, lda, B, ldb,
				    beta, C, ldc, H.kmax);
#endif
} // fgemm2

} // FFLAS

namespace FFLAS { namespace Protected{
	template <typename FloatElement, class Field>
	inline void  fgemm_convert (const Field& F,
				    const FFLAS_TRANSPOSE ta,
				    const FFLAS_TRANSPOSE tb,
				    const size_t m, const size_t n, const size_t k,
				    const typename Field::Element alpha,
				    const typename Field::Element* A,const size_t lda,
				    const typename Field::Element* B,const size_t ldb,
				    const typename Field::Element beta,
				    typename Field::Element * C, const size_t ldc,
				    const Winograd2Helper<FieldCategories::FloatingPointConvertibleTag>& H)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);
		FFPACK::UnparametricField<FloatElement> G;
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
		fgemm2(G, ta, tb, m, n, k, alphad,
		       Ad, ka, Bd, nb, betad, Cd, n,
		       Winograd2Helper<FieldCategories::FloatingPointTag>(H));
		// Conversion double = >  GFq
		finit(F, m, n, Cd, n, C, ldc);

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

	} // Protected
} // FFLAS


#endif // __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL
