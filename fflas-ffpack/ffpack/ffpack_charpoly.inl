/* ffpack/ffpack_charpoly.inl
 * Copyright (C) 2005 Clement Pernet
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
#ifndef __FFLASFFPACK_charpoly_INL
#define __FFLASFFPACK_charpoly_INL

#include "fflas-ffpack/utils/fflas_randommatrix.h"

namespace FFPACK {


    // 	template <class FloatElement, class Field, class Polynomial>
    // 	std::list<typename Polynomial<Element> >&
    // 	CharPoly_convert (const Field& F, std::list<typename Polynomial<Element> >& charp, const size_t N,
    // 		  typename Field::Element_ptr A, const size_t lda,
    // 		  const FFPACK_CHARPOLY_TAG CharpTag)
    // 	{
    // 		Givaro::ModularBalanced<FloatElement> G((FloatElement) F.cardinality());
    // 		FloatElement* Af = FFLAS::fflas_new<FloatElement>(N*N);
    // 		typename std::list< Polynomial<FloatElement> > charp_float;
    // 		fconvert(F, M, N, Af, N, A, lda);
    // // also convert the polynomial
    // 		CharPoly (G, charp_float, N, Af, N, CharpTag);

    // 		finit(F, ma, Yf, 1, Y, incY);
    // 		fflas_delete (Af);
    // 		return charp;
    // 	}
    template <class PolRing>
    std::list<typename PolRing::Element>&
    CharPoly (const PolRing& R, std::list<typename PolRing::Element>& charp, const size_t N,
              typename PolRing::Domain_t::Element_ptr A, const size_t lda,
              typename PolRing::Domain_t::RandIter& G,const FFPACK_CHARPOLY_TAG CharpTag,
              const size_t degree)
    {
        // if (Protected::AreEqual<PolRing::Domain_t, Givaro::Modular<double> >::value ||
        //     Protected::AreEqual<PolRing::Domain_t, Givaro::ModularBalanced<double> >::value){
        // 	if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
        // 		return CharPoly_convert <float,PolRing::Domain_t> (F, charp, N, A, lda, CharpTag);
        // }
        typedef typename PolRing::Domain_t Field;
        const Field& F = R.getdomain();

        FFPACK_CHARPOLY_TAG tag = CharpTag;
        if (tag == FfpackAuto){
            if (N < __FFLASFFPACK_CHARPOLY_Danilevskii_LUKrylov_THRESHOLD)
                tag = FfpackDanilevski;
            else if (N < __FFLASFFPACK_CHARPOLY_LUKrylov_ArithProg_THRESHOLD)
                tag = FfpackLUK;
            else
                tag = FfpackArithProgKrylovPrecond;
        }
        switch (tag){
            case FfpackDanilevski: return Danilevski (F, charp, N, A, lda);
            case FfpackLUK:
            {
                typename Field::Element_ptr X = FFLAS::fflas_new (F, N, N+1);
                Protected::LUKrylov (F, charp, N, A, lda, X, N, G);
                FFLAS::fflas_delete (X);
                return charp;
            }
            case FfpackArithProgKrylovPrecond:
            {
                size_t attempts=0;
                bool cont;
                
                Givaro::Integer p = F.characteristic();
                if (p < (uint64_t)N)	// Heuristic condition (the pessimistic theoretical one being p<2n^2).
                    return CharPoly(R, charp, N, A, lda, G, FfpackLUK);
                do{
                    typename Field::Element_ptr B = nullptr;
                    cont=false;
                    try {
                            // Preconditionning by a random block Krylov matrix.
                            // Some invariant factors may be discovered in the process and are stored in charp.
                        size_t ldb, Nb;
                        Protected::RandomKrylovPrecond (R, charp, N, A, lda, Nb, B, ldb, G, degree);
                            // Calling the main algorithm on the preconditionned part
                        Protected::ArithProg (R, charp, Nb, B, ldb, degree);
                        FFLAS::fflas_delete(B);
                    }
                    catch (CharpolyFailed){
                        if (B != nullptr)
                            FFLAS::fflas_delete(B);
                        charp.clear();
                        if (++attempts < 2)
                            cont = true;
                        else
                            return CharPoly (R, charp, N, A, lda, G, FfpackLUK);
                    }
                } while (cont);
                return charp;
            }
            case FfpackArithProg:
            {
                return Protected::ArithProg (R, charp, N, A, lda, 1);
            }

        case FfpackKG: return Protected::KellerGehrig (F, charp, N, A, lda);
        case FfpackKGFast:
                       {
                           size_t mc, mb, j;
                           if (Protected::KGFast (F, charp, N, A, lda, &mc, &mb, &j)){
                               std::cerr<<"NON GENERIC MATRIX PROVIDED TO KELLER-GEHRIG-FAST"<<std::endl;
                           }
                           return charp;
                       }
        case FfpackHybrid:
                       {
                           typename Field::Element_ptr X = FFLAS::fflas_new (F, N, N+1);
                           Protected::LUKrylov_KGFast (F, charp, N, A, lda, X, N);
                           FFLAS::fflas_delete (X);
                           return charp;
                       }
        case FfpackKGFastG:	return Protected::KGFast_generalized (F, charp, N, A, lda);
        default:
                            {
                                typename Field::Element_ptr X = FFLAS::fflas_new (F, N, N+1);
                                Protected::LUKrylov (F, charp, N, A, lda, X, N, G);
                                FFLAS::fflas_delete (X);
                                return charp;
                            }
        }
    }


    template <class PolRing>
    typename PolRing::Element&
    CharPoly (const PolRing& R, typename PolRing::Element& charp, const size_t N,
              typename PolRing::Domain_t::Element_ptr A, const size_t lda,
              typename PolRing::Domain_t::RandIter& G, const FFPACK_CHARPOLY_TAG CharpTag,
              const size_t degree){

        typedef typename PolRing::Domain_t Field;
        typedef typename PolRing::Element Polynomial;

        if (N==0)
            return R.assign (charp, R.one);

        Checker_charpoly<Field,Polynomial> checker(R.getdomain(),N,A,lda);

        std::list<Polynomial> factor_list;
        CharPoly (R, factor_list, N, A, lda, G, CharpTag, degree);
        typename std::list<Polynomial>::const_iterator it;
        it = factor_list.begin();

        R.init(charp, Givaro::Degree(N));
        R.assign (charp, *(it++));
        while( it!=factor_list.end() ){
            R.mulin(charp, *it);
            it++;
        }

        checker.check(charp);

        return charp;
    }


    namespace Protected {
        template <class Field, class Polynomial, class RandIter>
        std::list<Polynomial>&
        LUKrylov (const Field& F, std::list<Polynomial>& charp, const size_t N,
                  typename Field::Element_ptr A, const size_t lda,
                  typename Field::Element_ptr X, const size_t ldx, RandIter& G)
        {
            typedef typename Field::Element elt;
            elt* Ai, *Xi, *X2=X;
            size_t Ncurr=N;
            charp.clear();
            size_t nbfac = 0;

            while (Ncurr > 0){

                size_t *P = FFLAS::fflas_new<size_t>(Ncurr);
                typename Field::Element_ptr v = FFLAS::fflas_new (F,Ncurr,1);
                Polynomial minP;//=new Polynomial();
                //Hybrid_KGF_LUK_MinPoly (F, minP, (size_t)Ncurr, A, lda, X2, ldx, P);

                FFPACK::NonZeroRandomMatrix (F, 1, Ncurr, v, Ncurr, G);

                MatVecMinPoly (F, minP, Ncurr, A, lda, v, 1, X2, ldx, P);

                FFLAS::fflas_delete (v);

                size_t k = minP.size()-1; // degre of minpoly
                if ((k==1) && F.isZero (minP[0])){ // minpoly is X
                    if (FFLAS::fiszero(F,Ncurr, Ncurr, A, lda)){
                        // A is 0, CharPoly=X^n
                        minP.resize(Ncurr+1);
                        minP[1] = F.zero;
                        minP[Ncurr] = F.one;
                        k = Ncurr;
                    }
                }
                nbfac++;
                charp.push_front (minP);
                if (k==Ncurr){
                    FFLAS::fflas_delete( P);
                    return charp;
                }
                size_t Nrest = Ncurr-k;
                elt * X21 = X2 + k*ldx;
                elt * X22 = X21 + k;
                // Compute the n-k last rows of A' = PA^tP^t in X2_
                // A = A . P^t
                applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                        Ncurr, 0, k, A, lda, P);
                // Copy X2_ = (A'_2)^t
                for (Xi = X21, Ai = A+k; Xi != X21 + Nrest*ldx; Ai++, Xi+=ldx)
                    FFLAS::fassign(F, Ncurr, Ai, lda, Xi, 1);
                // A = A . P : Undo the permutation on A
                applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                        Ncurr, 0, k, A, lda, P);
                // X2_ = X2_ . P^t (=  (P A^t P^t)2_)
                applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                        Nrest, 0, k, X21, ldx, P);
                FFLAS::fflas_delete( P );
                // X21 = X21 . S1^-1
                ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper,
                      FFLAS::FflasNoTrans, FFLAS::FflasUnit, Nrest, k,
                      F.one, X2, ldx, X21, ldx);
                // Creation of the matrix A2 for recurise call
                for (Xi = X22, Ai = A;
                     Xi != X22 + Nrest*ldx;
                     Xi += (ldx-Nrest), Ai += (lda-Nrest))
                    for (size_t jj=0; jj<Nrest; ++jj)
                        *(Ai++) = *(Xi++);
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, k, F.mOne,
                       X21, ldx, X2+k, ldx, F.one, A, lda);
                X2 = X22;
                Ncurr = Nrest;
            }
            return charp;
        }


        template <class Field, class Polynomial>
        std::list<Polynomial>&
        LUKrylov_KGFast (const Field& F, std::list<Polynomial>& charp, const size_t N,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr X, const size_t ldx)
        {

            size_t kg_mc, kg_mb, kg_j;

            if (!KGFast (F, charp, N, A, lda, &kg_mc, &kg_mb, &kg_j))
                return charp;
            else{// Matrix A is not generic
                Polynomial *minP = new Polynomial();
                typename Field::ConstElement_ptr Ai;
                typename Field::Element_ptr A2i, Xi;
                size_t *P = FFLAS::fflas_new<size_t>(N);

                FFPACK::Protected::Hybrid_KGF_LUK_MinPoly (F, *minP, N, A, lda, X, ldx, P, FfpackKGF, kg_mc, kg_mb, kg_j);
                size_t k = minP->size()-1; // degre of minpoly
                if ((k==1) && F.isZero ((*minP)[0])){ // minpoly is X
                    Ai = A;
                    int j = int(N*N);
                    while (j-- && F.isZero(*(Ai++))) ;
                    if (!j){ // A is 0, CharPoly=X^n
                        minP->resize((size_t)N+1);
                        (*minP)[1] = F.zero;
                        (*minP)[N] = F.one;
                        k=N;
                    }
                }

                if (k==N){
                    charp.clear();
                    charp.push_front(*minP); // CharPoly = MinPoly
                    FFLAS::fflas_delete( P);
                    return charp;
                }

                size_t Nrest = N-k;
                typename Field::Element_ptr X21 = X + k*ldx;
                typename Field::Element_ptr X22 = X21 + k;

                // Creates the matrix A
                //size_t lambda = std::max(0,N - kg_mc*(kg_j+1) - kg_mb);  // uint >= 0 !!!
                size_t lambda =   kg_mc*(kg_j+1) + kg_mb;
                if (lambda > N)
                    lambda = 0 ;
                else
                    lambda = N - lambda ;

                size_t imax = kg_mc+kg_mb;
                // First Id
                for (size_t j = 0; j < lambda; ++j){
                    for (size_t i=0; i<imax; ++i)
                        F.assign (*(A+j+i*lda), F.zero);
                    F.assign (*(A+j+imax*lda), F.one);
                    for (size_t i=imax+1; i<N; ++i)
                        F.assign (*(A+j+i*lda), F.zero);
                    ++imax;
                }
                // Column block B
                for (typename Field::Element_ptr Aj=A; Aj<A+N*lda; Aj+=lda)
                    FFLAS::fassign (F, kg_mb, Aj+N-kg_mc-kg_mb, 1, Aj+lambda, 1);

                // Second Id block
                imax = N- kg_j*kg_mc;
                for (size_t j = 0; j< kg_j*kg_mc; ++j){
                    for (size_t i = 0; i<imax; ++i)
                        F.assign (*(A+lambda+kg_mb+j+i*lda), F.zero);
                    F.assign (*(A+lambda+kg_mb+j+imax*lda), F.one);
                    for (size_t i = imax+1; i<N; ++i)
                        F.assign (*(A+lambda+kg_mb+j+i*lda), F.zero);
                    ++imax;
                }

                // Compute the n-k last rows of A' = PA^tP^t in X2_

                // A = P . A
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                        N, 0,(int) k,
                        /*const_cast<typename Field::Element_ptr &>*/(A), lda, P);

                // Copy X2_ = (A'2_)
                for (Xi = X21, Ai = A+k*lda; Xi != X21 + Nrest*ldx; Ai+=lda-N, Xi+=ldx-N){
                    for (size_t jj=0; jj<N; ++jj){
                        F.assign(*(Xi++), *(Ai++));
                    }
                }

                // A = P^t . A : Undo the permutation on A
                applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                        N, 0,(int) k,
                        /*const_cast<typename Field::Element_ptr &>*/(A), lda, P);

                // X2_ = X2_ . P^t (=  (P A P^t)2_)
                applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                        Nrest, 0,(int) k, X21, ldx, P);

                // X21 = X21 . S1^-1
                ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Nrest, k,
                      F.one, X, ldx, X21, ldx);

                // Creation of the matrix A2 for recurise call
                typename Field::Element_ptr A2 = FFLAS::fflas_new (F, Nrest, Nrest);

                for (Xi = X22, A2i = A2;
                     Xi != X22 + Nrest*ldx;
                     Xi += (ldx-Nrest)){
                    for (size_t jj=0; jj<Nrest; ++jj){
                        *(A2i++) = *(Xi++);
                    }
                }
                fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, k, F.mOne,
                       X21, ldx, X+k, ldx, F.one, A2, Nrest);

                // Recursive call on X22
                LUKrylov_KGFast (F, charp, Nrest, A2, Nrest, X22, ldx);
                charp.push_front (*minP);
                FFLAS::fflas_delete( P);
                FFLAS::fflas_delete (A2);
                return charp;
            }
        }

    } // Protected

} // FFPACK
#include "ffpack_charpoly_mp.inl"
#endif // __FFLASFFPACK_charpoly_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
