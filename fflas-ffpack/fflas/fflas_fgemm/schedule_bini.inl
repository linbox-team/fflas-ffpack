/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 the LinBox group
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

/** @file fflas/fflas_fgemm/bini.inl
 * @ingroup MMalgos
 * @brief Bini implementation
 */

#ifndef __FFLASFFPACK_fgemm_bini_INL
#define __FFLASFFPACK_fgemm_bini_INL

namespace FFLAS { namespace BLAS3 {

	template < class Field >
	inline void Bini (const Field& F,
			  const FFLAS_TRANSPOSE ta,
			  const FFLAS_TRANSPOSE tb,
			  const size_t mr, const size_t nr, const size_t kr,
			  const typename Field::Element alpha,
			  const typename Field::Element* A,const size_t lda,
			  const typename Field::Element* B,const size_t ldb,
			  const typename Field::Element  beta,
			  typename Field::Element * C, const size_t ldc,
			  const size_t kmax, const size_t w, const FFLAS_BASE base,
			  const size_t rec_level)
	{

		FFLASFFPACK_check(F.isZero(beta));
		FFLASFFPACK_check(rec_level>0);

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

		namespace Protected {
#error "cacabouda"
			// C = a*A + B
			void add(const size_t m, const size_t n,
				 double a,
				 const double *A, const size_t lda,
				 const double *B, const size_t ldb,
				 double *C, const size_t ldc)
			{
				const double *Ai = A,*Bi = B;
				double *Ci       = C;
				for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
					for (size_t j = 0 ; j < n ; ++j)
						Ci[j] = a * Ai[j] + Bi[j];
			}

			// C = C-(A+B)
			void subadd(const size_t m, const size_t n,
				    const double *A, const size_t lda,
				    const double *B, const size_t ldb,
				    double *C, const size_t ldc)
			{
				const double *Ai = A,*Bi = B;
				double *Ci       = C;
				for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
					for (size_t j = 0 ; j < n ; ++j) {
						Ci[j] = Ci[j] - Ai[j] - Bi[j] ;
					}

			}

			// C = C+A-B
			void addsub(const size_t m, const size_t n,
				    const double *A, const size_t lda,
				    const double *B, const size_t ldb,
				    double *C, const size_t ldc)
			{
				const double *Ai = A,*Bi = B;
				double *Ci       = C;
				for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
					for (size_t j = 0 ; j < n ; ++j) {
						Ci[j] = Ci[j] + Ai[j] - Bi[j] ;
					}

			}


			// C = (C+B)/e
			template<class Field>
			void addscalinf(const Field & F, const size_t m, const size_t n,
					const double *B, const size_t ldb,
					double e,
					double *C, const size_t ldc)
			{
				const double * Bi = B;
				double * Ci = C;
				for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
					for (size_t j = 0 ; j < n ; ++j)
						F.init( Ci[j], (Ci[j]+Bi[j])/e );

			}

			// C = (C-B)/e
			template<class Field>
			void subscalinf(const Field & F, const size_t m, const size_t n,
					const double *B, const size_t ldb,
					double e,
					double *C, const size_t ldc)
			{
				const double * Bi = B;
				double * Ci = C;
				for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
					for (size_t j = 0 ; j < n ; ++j)
						F.init( Ci[j], (Ci[j]-Bi[j])/e );

			}


			template<class Field>
			double * gemm_fflas(const Field & F,
					    const size_t m, const size_t n, const size_t k,
					    const double *A, size_t lda,
					    const double *B, size_t ldb,
					    double * C, size_t ldc,
					    int rec =  0)
			{
				FFLAS::DoubleDomain R;
				FFLAS::fgemm(R,
					     FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
					     m,n,k,
					     1,
					     A,lda, B,ldb,
					     0,
					     C, ldc);
				// Field F2 (F.characteristic()*F.characteristic());
				// FFLAS::finit(F2,m,n,C,ldc);
				return C;
			}
		} // Protected


			namespace Protected {
				namespace Rec {
					// Field must be Modular<double>
					template<class Field>
					double *
					gemm_bini(const Field & F
						  , const size_t m
						  , const size_t n
						  , const size_t k
						  , const double *A , const size_t lda
						  , const double *B , const size_t ldb
						  , double *C , const size_t ldc
						  , int rec
						 )
					{
						FFPACK::UnparametricField<double>   NoField;
						const double p = (double)F.characteristic();
						size_t M = (n>m)?std::min(k,m):std::min(k,n);
						// std::cout << rec << ',' <<  M  << std::endl;

						if ( M <= TRE  || rec <= 0) {
							// std::cout << "ffw" << std::endl;
							return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
							// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
						}

						assert(k/2*2==k); // k divisible par 2
						assert(n/2*2==n); // n divisible par 2
						assert(m/3*3==m); // m divisible par 3


						size_t n2 = n/2;
						size_t k2 = k/2;
						size_t m3 = m/3;

						const double epsilon = p ;
						// const double epsilon = 1/(double)(1<<26);
						// std::cout << "€ = " << epsilon << std::endl;

						// sub matrices in A
						const double * A11 = A;
						const double * A12 = A   +k2;
						const double * A21 = A   +lda*m3;
						const double * A22 = A21 +k2;
						const double * A31 = A21 +lda*m3;
						const double * A32 = A31 +k2;

						// sub matrices in C
						double * C11 = C;
						double * C12 = C   +n2;
						double * C21 = C   +ldc*m3;
						double * C22 = C21 +n2;
						double * C31 = C21 +ldc*m3;
						double * C32 = C31 +n2;

						// sub matrices in B
						const double * B11 = B;
						const double * B12 = B   +n2;
						const double * B21 = B   +ldb*k2;
						const double * B22 = B21 +n2;

						FFLAS::fzero(F,m,n,C,ldc);

						/*
						 * Algo :
						 * S1  := A11  +A22;
						 * S4  := e*A12+A22;
						 * S5  := A11  +e*A12;
						 * S6  := A21  +A32;
						 * S9  := A21  +e*A31;
						 * S10 := e*A31+A32;
						 *
						 * T1  := e*B11 +B22;
						 * T2  := B21   +B22;
						 * T4  := -e*B11+B21;
						 * T5  := e*B12 +B22;
						 * T6  := B11   +e*B22;
						 * T7  := B11   +B12;
						 * T9  := B12   -e*B22;
						 * T10 := B11   +e*B21;
						 *
						 * P1 := S1 *T1;
						 * P2 := A22*T2;
						 * P3 := A11*B22;
						 * P4 := S4 *T4;
						 * P5 := S5 *T5;
						 * P6 := S6 *T6;
						 * P7 := A21*T7;
						 * P8 := A32*B11;
						 * P9 := S9 *T9;
						 * P10:= S10*T10;
						 *
						 * C11 := (P1-P2-P3+P4)/e;
						 * C12 := (P3-P5)/(-e) ;
						 * C21 := P4+P6-P10 ;
						 * C22 := P1-P5+P9;
						 * C31 := (-P8+P10)/e;
						 * C32 := (P6-P7-P8+P9)/e;
						 *
						 */

						double * eA12 = new double [m3*k2];
						double * eA31 = eA12 ;

						double * S1 = new double[m3*k2] ;
						double * T1 = new double[n2*k2] ; // ou aire
						double * P1 = new double[n2*m3] ; // ou aire
						double * P2 = new double[n2*m3] ; // ou aire
						// double * C11t = new double[n2*m3] ;
						// S1  := A11  +A22;
						// XXX add(m3,k2,A11,lda,A22,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,S1,k2);
						// T1  := e*B11 +B22;
						add(k2,n2,epsilon,B11,ldb,B22,ldb,T1,n2);
						// P1 := S1 *T1; (dans C22)
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,C22,ldc,rec-1);
						// S4  := e*A12+A22;
						// XXX scal(m3,k2,epsilon,eA12,k2,A12,lda) ;
						FFLAS::fscal(NoField,m3,k2,epsilon,eA12,k2,A12,lda) ;
						// XXX add(m3,k2,eA12,k2,A22,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,eA12,k2,A22,lda,S1,k2);
						// T4  := -e*B11+B21;
						add(k2,n2,-epsilon,B11,ldb,B21,ldb,T1,n2);
						// P4 := S4 *T4; (dans C21)
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,C21,ldc,rec-1);
						// C11 = P1+P4
						// XXX add(m3,n2,C21,ldc,C22,ldc,C11,ldc);
						FFLAS::fadd(NoField,m3,n2,C21,ldc,C22,ldc,C11,ldc);
						// T2  := B21  +B22;
						// XXX add(k2,n2,B21,ldb,B22,ldb,T1,n2);
						FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,T1,n2);
						// P2 := A22*T2;
						gemm_bini(F,m3,n2,k2,A22,lda,T1,n2,P1,n2,rec-1);
						// P3 := A11*B22; (dans C12)
						gemm_bini(F,m3,n2,k2,A11,lda,B22,ldb,C12,ldc,rec-1);
						// C11 -= (P2+P3)
						subadd(m3,n2,P1,n2,C12,ldc,C11,ldc);
						// S5  := A11  +e*A12;
						// XXX add(m3,k2,eA12,k2,A11,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,eA12,k2,A11,lda,S1,k2);
						// T5  := e*B12 +B22;
						add(k2,n2,epsilon,B12,ldb,B22,ldb,T1,n2);
						// P5 := S5 *T5;
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,P2,n2,rec-1);
						// C12 -= P5
						// XXX subscalinf(p,m3,n2,P2,n2,-epsilon,C12,ldc);
						subscalinf(NoField,m3,n2,P2,n2,-epsilon,C12,ldc);
						// S6  := A21  +A32;
						// XXX add(m3,k2,A21,lda,A32,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,S1,k2);
						// T6  := B11   +e*B22;
						add(k2,n2,epsilon,B22,ldb,B11,ldb,T1,n2);
						// P6 := S6 *T6; dans C32
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,C32,ldc,rec-1);
						// C21+= P6
						// XXX addin(m3,n2,C32,ldc,C21,ldc);
						FFLAS::faddin(NoField,m3,n2,C32,ldc,C21,ldc);
						// T7  := B11  +B12;
						// XXX add(k2,n2,B11,ldb,B12,ldb,T1,n2);
						FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,T1,n2);
						// P7 := A21*T7; !signe
						gemm_bini(F,m3,n2,k2,A21,lda,T1,n2,P1,n2,rec-1);
						// P8 := A32*B11; dans C31 !signe
						gemm_bini(F,m3,n2,k2,A32,lda,B11,ldb,C31,ldc,rec-1);
						// C32 -= P8+P7
						subadd(m3,n2,P1,n2,C31,ldc,C32,ldc);
						// S9  := A21  +e*A31;
						// XXX scal(m3,k2,epsilon,eA31,k2,A31,lda);
						FFLAS::fscal(NoField,m3,k2,epsilon,eA31,k2,A31,lda);
						// XXX add(m3,k2,eA31,k2,A21,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,eA31,k2,A21,lda,S1,k2);
						// T9  := B12   -e*B22;
						add(k2,n2,-epsilon,B22,ldb,B12,ldb,T1,n2);
						// P9 := S9 *T9;
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1);
						// C32= (C32+P9)/p
						// XXX addscalinf(p,m3,n2,P1,n2,epsilon,C32,ldc);
						addscalinf(NoField,m3,n2,P1,n2,epsilon,C32,ldc);
						// C22+= P9-P5
						addsub(m3,n2,P1,n2,P2,n2,C22,ldc);
						// S10 := e*A31+A32;
						// XXX add(m3,k2,eA31,k2,A32,lda,S1,k2);
						FFLAS::fadd(NoField,m3,k2,eA31,k2,A32,lda,S1,k2);
						// T10 := B11   +e*B21;
						add(k2,n2,epsilon,B21,ldb,B11,ldb,T1,n2);
						// P10:= S10*T10;
						gemm_bini(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1);
						// C21-= P10
						// XXX subin(m3,n2,P1,n2,C21,ldc);
						FFLAS::fsubin(NoField,m3,n2,P1,n2,C21,ldc);
						// C31= (C31-P10)/(-epsilon)
						// XXX subscalinf(p,m3,n2,P1,n2,-epsilon,C31,ldc);
						subscalinf(NoField,m3,n2,P1,n2,-epsilon,C31,ldc);
						// C11 := (P1+P-P3+P4)/e;
						// XXX scalinf(p,m3,n2,epsilon,C11,ldc);
						// XXX FFLAS::fscalin(F,m3,n2,(double)1/epsilon,C11,ldc);
						FFLAS::fscalin(NoField,m3,n2,(double)1/epsilon,C11,ldc);

						// FFLAS::finit(F2,m,n,C,ldc);
						// Field F2 (F.characteristic()*F.characteristic());
						// FFLAS::finit(F2,m,n,C,ldc);

						delete[] P1;
						delete[] P2;
						delete[] S1;
						delete[] T1;
						delete[] eA12 ;

						return C;

					}

				} // Rec

				template<class Field>
				typename Field::Element *
				gemm_bini(const Field &F
					  , const size_t m
					  , const size_t n
					  , const size_t k
					  , const typename Field::Element *A
					  , const size_t lda
					  , const typename Field::Element *B
					  , const size_t ldb
					  , typename Field::Element *C
					  , const size_t ldc
					  , int rec
					 )
				{
					std::cout << "threshold = " << TRE << ", rec level =" << rec << std::endl;

					assert(k/2*2==k); // k divisible par 2
					assert(n/2*2==n); // n divisible par 2
					assert(m/3*3==m); // m divisible par 3

					// e-formule
					Rec::gemm_bini(F,m,n,k,A,lda,B,ldb,C,ldc,rec);

					// vire les e.
					FFLAS::finit(F,m,n,C,ldc);

					return C;

				}

			}


		} // Winograd

	} // BLAS3


} // FFLAS

#endif // __FFLASFFPACK_fgemm_bini_INL

