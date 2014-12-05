/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* ffpack/ffpack_charpoly_kgfast.inl
 * Copyright (C) 2004 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_ffpack_charpoly_kgfast_INL
#define __FFLASFFPACK_ffpack_charpoly_kgfast_INL

namespace FFPACK {
	namespace Protected {
		//---------------------------------------------------------------------
		// CharPoly: Compute the characteristic polynomial of A using
		// Keller-Gehrig's fast algorithm. A must be generic.
		//---------------------------------------------------------------------

		template <class Field, class Polynomial>
		int
		KGFast ( const Field& F, std::list<Polynomial>& charp,
			 const size_t N,
			 typename Field::Element_ptr A, const size_t lda,
			 size_t * kg_mc, size_t* kg_mb, size_t* kg_j )
		{

			//std::cerr<<"Dans KGFast"<<std::endl;
			size_t mc=N>>1; // Matrix A is transformed into a mc_Frobenius form
			size_t mb=N-mc;
			typename Field::Element_ptr C, B;


			while ( mc > 0 ) {
				// size_t r;
#if 0
				std::cerr<<"Boucle1: mc,mb,N="<<mc<<" "<<mb<<" "<<N<<std::endl;
				write_field( F, std::cerr, A, N, N, lda );
#endif
				size_t j=0;
				C = A + (N-mc);
				//std::cerr<<std::endl<<"mc="<<mc<<":";
				while ( (j+1)*mc < N ) {
					mb = std::min ( mb, N-(j+1)*mc );
#if 0
					std::cerr<<"Boucle2: j,mb="<<j<<" "<<mb<<std::endl;
					write_field( F, std::cerr, A, N, N, lda );
#endif
					B = A + (N-mc-mb);

					// B1 <- C1^-1.B1
					typename Field::Element_ptr LUP = FFLAS::fflas_new (F, mc, mc);
					// for (size_t i=0; i<mc; ++i)
						// FFLAS::fassign( F, mc, C+i*lda, 1, LUP+i*mc, 1);
					FFLAS::fassign(F,mc,mc,C,lda,LUP,mc);
					size_t * P = FFLAS::fflas_new<size_t>(mc);
					size_t * Q = FFLAS::fflas_new<size_t>(mc);

					if ( (LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mc, mc, LUP, mc, P, Q)) < mc ){
						* kg_mc = mc;
						* kg_mb = mb;
						* kg_j = j;
						FFLAS::fflas_delete( P);
						FFLAS::fflas_delete( Q);
						FFLAS::fflas_delete (LUP);
						return -1;

					}
#if 0
					std::cerr<<"LUP="<<std::endl;
					write_field( F, std::cerr, LUP, mc, mc, mc );
#endif
					//std::cerr<<" "<<r;
					ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
					      mc, mb, F.one, LUP, mc , B, lda);
					ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      mc, mb, F.one, LUP, mc , B, lda);
					FFLAS::fflas_delete (LUP);
					applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, mb, 0, (int)mc, B, lda, P );

					FFLAS::fflas_delete( P);
					FFLAS::fflas_delete( Q);
#if 0
					std::cerr<<"Apres B1<-C1^-1"<<std::endl;
					write_field( F, std::cerr, A, N, N, lda );
#endif

					// B2 <- B2 - C2.B1
					fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-mc, mb, mc,
					      F.mOne, C+mc*lda, lda, B, lda,
					      F.one, B+mc*lda, lda);

#if 0
					std::cerr<<"Apres B2<-B2-C2.B1"<<std::endl;
					write_field( F, std::cerr, A, N, N, lda );
#endif

					// Shifting B: B1;B2 -> B2;B1
					typename Field::Element_ptr tmp = FFLAS::fflas_new (F, mc, mb);
					// for (size_t i=0; i<mc; ++i)
						// FFLAS::fassign( F, mb, B+i*lda, 1, tmp+i*mb, 1);
					FFLAS::fassign(F,mc,mb,B,lda,tmp,mb);
					// for (size_t i=mc; i<N; ++i)
						// FFLAS::fassign( F, mb, B+i*lda, 1, B+(i-mc)*lda, 1);
					FFLAS::fassign(F,N-mc,mb,B+mc*lda,lda,B,lda);
					// for (size_t i=0; i<mc; ++i)
						// FFLAS::fassign( F, mb, tmp+i*mb, 1, B+(i+N-mc)*lda, 1);
					FFLAS::fassign(F,mc,mb,tmp,mb,B+(N-mc)*lda,lda);
					FFLAS::fflas_delete (tmp);
#if 0
					std::cerr<<"Apres shift de B"<<std::endl;
					write_field( F, std::cerr, A, N, N, lda );
#endif


					// C3 <- B3.C1 + C3
					fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, (j+1)*mc, mc, mb,
					      F.one, B+(N-(j+1)*mc)*lda, lda, C+(N-(j+1)*mc-mb)*lda, lda,
					      F.one, C+(N-(j+1)*mc)*lda, lda);
#if 0
					std::cerr<<"C3 <- B3.C1 + C3: B3="<<std::endl;
					write_field( F, std::cerr, B+(N-(j+1)*mc)*lda, (j+1)*mc, mb, lda );
					std::cerr<<"C3 <- B3.C1 + C3: C1"<<std::endl;
					write_field( F, std::cerr,  C+(N-(j+1)*mc-mb)*lda, mb, mc, lda );
					std::cerr<<"C3 <- B3.C1 + C3: C3="<<std::endl;
					write_field( F, std::cerr, C+(N-(j+1)*mc)*lda, (j+1)*mc, mc, lda );
#endif

					int lambda = (int)(N - mb - (j+1)*mc );
					if ( int(mb) < lambda ){

#if 0
						std::cerr<<"mb<lambda"<<std::endl;
#endif
						typename Field::Element_ptr tmp2 = FFLAS::fflas_new (F, (size_t)lambda, mc);

						// tmp2 <- C1
						// for (int i=0; i<lambda; ++i)
							// FFLAS::fassign( F, mc, C+i*(int)lda, 1, tmp2+i*(int)mc, 1);
						FFLAS::fassign(F,(size_t)lambda,mc,C,lda,tmp2,mc);

						// C1' <- B1.C2
						fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mb, mc, mb,
						      F.one, B, lda, C+lambda*(int)lda, lda,
						      F.zero, C, lda);

						// tmp2 <- B2.C2 + tmp2
						fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, (size_t)lambda, mc, mb,
						      F.one, B+mb*lda, lda, C+lambda*(int)lda, lda,
						      F.one, tmp2, mc);

						// C2' <- tmp2
						// for (int i=0; i<lambda; ++i)
							// FFLAS::fassign( F, mc, tmp2+(size_t)i*mc, 1, C+(mb+(size_t)i)*lda, 1);
						FFLAS::fassign(F,(size_t)lambda,mc,tmp2,mc,C+mb*lda,lda);
						FFLAS::fflas_delete (tmp2);
					}
					else if ( lambda > 0 ){
#if 0
						std::cerr<<"lambda>0"<<std::endl;
#endif

						typename Field::Element_ptr tmp2 = FFLAS::fflas_new (F, mb, mc);
						// C1 <- B2.C2 + C1
						fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, (size_t)lambda, mc, mb,
						      F.one, B+mb*lda, lda, C+lambda*(int)lda, lda,
						      F.one, C, lda);

						// tmp2 <-B1.C2
						fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mb, mc, mb,
						      F.one, B, lda, C+lambda*(int)lda, lda,
						      F.zero, tmp2, mc);

						// C2' <- C1
						// for (int i=0; i<lambda; ++i)
							// FFLAS::fassign( F, mc, C+i*(int)lda, 1, C+(mb+(size_t)i)*lda, 1);
						FFLAS::fassign(F,(size_t)lambda,mc,C,lda,C+mb*lda,lda);

						// C1' <- tmp2
						// for (size_t i=0; i<mb; ++i)
							// FFLAS::fassign( F, mc, tmp2+i*mc, 1, C+i*lda, 1);
						FFLAS::fassign(F,mb,mc,tmp2,mc,C,lda);
						FFLAS::fflas_delete (tmp2);
					}
					else{
#if 0
						std::cerr<<"lambda<0"<<std::endl;
#endif
						mb = N - (j+1)*mc;
						typename Field::Element_ptr tmp2 = FFLAS::fflas_new (F, mb, mc);

						// tmp2 <-B1.C1
						fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mb, mc, mb,
						      F.one, B, lda, C, lda,
						      F.zero, tmp2, mc);

						// C1' <- tmp2
						// for (size_t i=0; i<mb; ++i)
							// FFLAS::fassign( F, mc, tmp2+i*mc, 1, C+i*lda, 1);
						FFLAS::fassign(F,mb,mc,tmp2,mc,C,lda);
						FFLAS::fflas_delete (tmp2);
					}

					j++;
				}
				mb = mc;
				mc>>=1;
				mb -= mc;
			}

			Polynomial *minP = new Polynomial();
			minP->resize(N+1);
			minP->operator[](N) = F.one;
			typename Polynomial::iterator it = minP->begin();
			for (size_t j=0; j<N; ++j, it++){
				F.neg(*it, *(A+N-1+j*lda));
			}
			charp.clear();
			charp.push_front(*minP);
			return 0;
		}

		template<class Field>
		void
		fgemv_kgf( const Field& F,  const size_t N,
			   typename Field::ConstElement_ptr A, const size_t lda,
			   typename Field::ConstElement_ptr X, const size_t incX,
			   typename Field::Element_ptr Y, const size_t incY,
			   const size_t kg_mc, const size_t kg_mb, const size_t kg_j )
		{

			size_t big_truc =kg_mb-kg_mc*(kg_j+1) ;
			size_t lambda = (N<big_truc)?(0):(N-big_truc);
			// Y1 <- X2
			FFLAS::fassign ( F, lambda, X+(kg_mb+kg_mc)*incX, incX, Y, incY );
			// Y2 <- X.B
			fgemv( F, FFLAS::FflasTrans, N, kg_mb, F.one, A+N-kg_mc-kg_mb, lda, X, incX, F.zero, Y+lambda*incY, incY );
			// Y3 <- X3
			FFLAS::fassign ( F, kg_j*kg_mc, X+(lambda+kg_mb+kg_mc)*incX, incX, Y+(lambda+kg_mb)*incY, incY );
			// Y4 <- X.C
			fgemv( F, FFLAS::FflasTrans, N, kg_mc, F.one, A+N-kg_mc, lda, X, incX, F.zero, Y+(N-kg_mc)*incY, incY );
		}

	} // Protected

} // FFPACK

#endif // __FFLASFFPACK_ffpack_charpoly_kgfast_INL
