/* ffpack/ffpack_minpoly.inl
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
#ifndef __FFLASFFPACK_ffpack_minpoly_INL
#define __FFLASFFPACK_ffpack_minpoly_INL
namespace FFPACK {

    template <class Field, class Polynomial>
    inline Polynomial&
    MinPoly (const Field& F, Polynomial& minP, const size_t N,
             typename Field::ConstElement_ptr A, const size_t lda){

        typename Field::RandIter G (F);
        return MinPoly(F, minP, N, A, lda, G);
    }

    template <class Field, class Polynomial, class RandIter>
    inline Polynomial&
    MinPoly (const Field& F, Polynomial& minP, const size_t N,
             typename Field::ConstElement_ptr A, const size_t lda,
             RandIter& G){

        if (N==0){
            minP.resize(1);
            F.assign(minP[0],F.one);
            return minP;
        }
        // Allocating a Krylov basis
        typename Field::Element_ptr v = FFLAS::fflas_new(F, 1, N);
        // Picking a non-zero random vector
        NonZeroRandomMatrix (F, 1, N, v, N, G);

        MatVecMinPoly (F, minP, N, A, lda, v, 1);

        FFLAS::fflas_delete(v);
        return minP;
    }

    template <class Field, class Polynomial>
    inline Polynomial&
    MatVecMinPoly (const Field& F, Polynomial& minP, const size_t N,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::ConstElement_ptr v, const size_t incv){

        // Construct the Krylov matrix K and eliminate it online
        typename Field::Element_ptr K = FFLAS::fflas_new(F, N+1, N);
        const size_t ldk = N;
        size_t * P = FFLAS::fflas_new<size_t>(N);

        typename Field::Element_ptr u = FFLAS::fflas_new(F, 1, N);
        FFLAS::fassign (F, N, v, incv, u, 1);

        Protected::MatVecMinPoly(F, minP, N, A, lda, u, 1, K, ldk, P);

        FFLAS::fflas_delete (u);
        FFLAS::fflas_delete (P);
        FFLAS::fflas_delete (K);

        return minP;
    }


    namespace Protected {

        template <class Field, class Polynomial>
        inline Polynomial&
        MatVecMinPoly (const Field& F, Polynomial& minP, const size_t N,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::Element_ptr v, const size_t incv,
                       typename Field::Element_ptr K, const size_t ldk,
                       size_t * P){

            FFLAS::fassign (F, N, v, incv, K, 1);

            // Construct the Krylov matrix K and eliminate it online
            // LUP factorization on K, construct K on the fly

            size_t k = Protected::LUdivine_construct (F, FFLAS::FflasUnit, N+1, N, A, lda, K, ldk, v, incv, P, true, FfpackDense);


            minP.resize(k+1);

            minP[k] = F.one;
            // If minpoly is X
            if (k==1 && F.isZero (*(K+ldk))){
                minP[0] = F.zero;
                return minP;
            }
            // Get the coefficients from the first linear dependency in K
            typename Field::Element_ptr Kk = K+k*ldk;
            ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasNonUnit, k, K, ldk, Kk, 1);
            for (size_t j=0; j<k; ++j)
                F.neg (minP[j],Kk[j]);
            return minP;
        }

        template <class Field, class Polynomial>
        Polynomial&
        Hybrid_KGF_LUK_MinPoly (const Field& F, Polynomial& minP, const size_t N
                                ,typename Field::ConstElement_ptr A, const size_t lda
                                ,typename Field::Element_ptr X, const size_t ldx
                                ,size_t* P
                                ,const FFPACK_MINPOLY_TAG MinTag// = FfpackDense
                                ,const size_t kg_mc// =0
                                ,const size_t kg_mb//=0
                                ,const size_t kg_j //=0
                               )
        {
            // nRow is the number of row in the krylov base already computed
            size_t j, k ;
            //size_t	nRow = 2;
            typename Polynomial::iterator it;
            typename Field::Element_ptr Xi, Ui;
            typename Field::RandIter g (F);
            bool KeepOn=true;
            typename Field::Element_ptr U = FFLAS::fflas_new (F, N);
            // Picking a non zero vector
            do{
                for (Ui=U, Xi = X; Ui<U+N; ++Ui, ++Xi){
                    g.random (*Ui);
                    *Xi = *Ui;
                    if (!F.isZero(*Ui))
                        KeepOn = false;
                }
            }while(KeepOn);

            // nRow = 1;
            // LUP factorization of the Krylov Base Matrix
            k = Protected::LUdivine_construct (F, FFLAS::FflasUnit, N+1, N, A, lda, X, ldx, U, 1, P, true,
                                               MinTag, kg_mc, kg_mb, kg_j);
            //FFLAS::fflas_delete( U);
            minP.resize(k+1);
            minP[k] = F.one;
            if ( (k==1) && F.isZero(*(X+ldx))){ // minpoly is X
                FFLAS::fflas_delete (U);
                for (size_t i=0; i<k; ++i)
                    minP[i] = F.zero;
                return minP;
            }
            // U contains the k first coefs of the minpoly
            //typename Field::Element_ptr m= FFLAS::fflas_new<elt>(k);
            FFLAS::fassign( F, k, X+k*ldx, 1, U, 1);
            ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasNonUnit, k, X, ldx, U, 1);
            it = minP.begin();
            for (j=0; j<k; ++j, it++){
                F.neg(*it, U[j]);
            }
            FFLAS::fflas_delete (U);
            return minP;
        }
    } // Protected

} // FFPACK
#endif // __FFLASFFPACK_ffpack_minpoly_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
