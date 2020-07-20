/* ffpack/ffpack_ludivine.inl
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

#ifndef __FFLASFFPACK_ffpack_ludivine_INL
#define __FFLASFFPACK_ffpack_ludivine_INL

#include "fflas-ffpack/fflas/fflas_bounds.inl"

namespace FFPACK {
    template<class Field>
    inline size_t
    LUdivine_gauss( const Field& F, const FFLAS::FFLAS_DIAG Diag,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda, size_t*P,
                    size_t *Q, const FFPACK::FFPACK_LU_TAG LuTag)
    {
        size_t MN = std::min(M,N);
        typename Field::Element_ptr Acurr = A;
        size_t r = 0;

        for (size_t k = 0; k < MN; ++k){
            size_t p = r;
            Acurr = A+k*lda+r;
            while ((p < N) && F.isZero (*(Acurr++)))
                p++;
            if (p < N){
                P[r] = p;
                if (r < k){
                    FFLAS::fassign (F, N-r, (A+k*lda+r),1, (A + r*(lda+1)), 1);
                    Acurr = A+r+k*lda;
                    for (size_t i=r; i<N; ++i)
                        F.assign(*(Acurr++),F.zero);
                }

                FFLAS::fswap (F, M, A+r, lda, A+p, lda);
                Q[r] = k;
                r++;
            }
            if (k+1<M){
                ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, r, A, lda, A+(k+1)*lda, 1);
                fgemv (F, FFLAS::FflasTrans, r, N-r, F.mOne, A+r, lda, A+(k+1)*lda, 1, F.one, A+(k+1)*lda+r, 1);
            }
            else
                break; // return r;
        }
        return r;
    }

    template<class Element>
    class callLUdivine_small;



    template<class Field>
    inline size_t
    LUdivine_small( const Field& F, const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda, size_t*P,
                    size_t *Q, const FFPACK::FFPACK_LU_TAG LuTag)
    {
        return callLUdivine_small <typename Field::Element> ()
        (F, Diag, trans, M, N, A, lda, P, Q, LuTag);
    }

    template<class Element>
    class callLUdivine_small {
    public:
        template <class Field>
        inline size_t
        operator()( const Field& F, const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda, size_t*P,
                    size_t *Q, const FFPACK::FFPACK_LU_TAG LuTag)
        {
            if ( !(M && N) ) return 0;
            typedef typename Field::Element elt;
            typedef typename Field::Element_ptr elt_ptr;
            elt_ptr Aini = A;
            elt_ptr Acurr;
            size_t rowp = 0;
            size_t R = 0;
            size_t k = 0;
            while ((rowp<M) && (k<N)){
                size_t colp;

                //Find non zero pivot
                colp = k;
                Acurr = Aini;
                while ((F.isZero(*Acurr)) || (F.isZero (F.reduce (*Acurr))))
                    if (++colp == N){
                        if (rowp==M-1)
                            break;
                        colp=k; ++rowp;
                        Acurr = Aini += lda;
                    }
                    else
                        ++Acurr;

                if ((rowp == M-1)&&(colp == N))
                    break;
                R++;
                P[k] = colp;
                Q[k] = rowp;

                // Permutation of the pivot column
                FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

                //Normalization
                elt invpiv;
                F.init (invpiv);
                F.reduce (*Aini);
                F.inv (invpiv,*Aini);

                for (size_t j=1; j<N-k; ++j)
                    if (!F.isZero(*(Aini+j)))
                        F.reduce (*(Aini+j));
                for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    if (!F.isZero(*(Aini+i)))
                        F.reduce (*(Aini+i));


                if (Diag == FFLAS::FflasUnit) {
                    // for (size_t j=1; j<N-k; ++j)
                    // if (!F.isZero(*(Aini+j)))
                    // F.mulin (*(Aini+j),invpiv);
                    FFLAS::fscalin(F,N-k-1,invpiv,Aini+1,1);
                }
                else {
                    // for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    // if (!F.isZero(*(Aini+i)))
                    // F.mulin (*(Aini+i),invpiv);
                    FFLAS::fscalin(F,M-rowp-1,invpiv,Aini+lda,lda);
                }

                //Elimination
                //Or equivalently, but without delayed ops :
                FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

                Aini += lda+1; ++rowp; ++k;
            }

            // Compression the U matrix
            size_t l;
            if (Diag == FFLAS::FflasNonUnit){
                Aini = A;
                l = N;
            }
            else {
                Aini = A+1;
                l=N-1;
            }
            for (size_t i=0; i<R; ++i, Aini += lda+1) {
                if (Q[i] > i){
                    FFLAS::fassign (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
                    for (size_t j=0; j<l-i; ++j)
                        F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
                }
            }
            return R;
        }
    };

    template<>
    class callLUdivine_small<double> {
    public:
        template <class Field>
        inline size_t
        operator()( const Field& F,
                    const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda, size_t*P,
                    size_t *Q, const FFPACK::FFPACK_LU_TAG LuTag)
        {

            if ( !(M && N) ) return 0;
            typedef typename Field::Element elt;
            elt * Aini = A;
            elt * Acurr;
            size_t rowp = 0;
            size_t R = 0;
            size_t k = 0;
            size_t delay =0;
            size_t kmax = FFLAS::Protected::DotProdBoundClassic (F, F.one) -1; // the max number of delayed operations
            while ((rowp<M) && (k<N)){
                size_t colp;

                //Find non zero pivot
                colp = k;
                Acurr = Aini;
                while ((F.isZero(*Acurr)) || (F.isZero (F.reduce (*Acurr))))
                    if (++colp == N){
                        if (rowp==M-1)
                            break;
                        colp=k; ++rowp;
                        Acurr = Aini += lda;
                    }
                    else
                        ++Acurr;

                if ((rowp == M-1)&&(colp == N))
                    break;
                R++;
                P[k] = colp;
                Q[k] = rowp;

                // Permutation of the pivot column
                FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

                //Normalization
                elt invpiv;
                F.reduce (*Aini);
                F.inv (invpiv,*Aini);

                for (size_t j=1; j<N-k; ++j)
                    if (!F.isZero(*(Aini+j)))
                        F.reduce (*(Aini+j));
                for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    if (!F.isZero(*(Aini+i)))
                        F.reduce(*(Aini+i));


                if (Diag == FFLAS::FflasUnit) {
                    // for (size_t j=1; j<N-k; ++j)
                    // if (!F.isZero(*(Aini+j)))
                    // F.mulin (*(Aini+j),invpiv);
                    FFLAS::fscalin(F,N-k-1,invpiv,Aini+1,1);
                }
                else {
                    // for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    // if (!F.isZero(*(Aini+i)))
                    // F.mulin (*(Aini+i),invpiv);
                    FFLAS::fscalin(F,M-rowp-1,invpiv,Aini+lda,lda);
                }

                if (delay++ >= kmax){ // Reduction has to be done
                    delay = 0;
                    FFLAS::freduce (F, M-rowp-1,N-k-1, Aini+lda+1, lda);
                    // for (size_t i=1; i<M-rowp; ++i)
                    // 	for (size_t j=1; j<N-k; ++j)
                    // 		F.init(	*(Aini+i*lda+j),*(Aini+i*lda+j));
                }
                //Elimination
                for (size_t i=1; i<M-rowp; ++i)
                    for (size_t j=1; j<N-k; ++j)
                        *(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
                //Or equivalently, but without delayed ops :
                //FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

                Aini += lda+1; ++rowp; ++k;
            }

            // Compression the U matrix
            size_t l;
            if (Diag == FFLAS::FflasNonUnit){
                Aini = A;
                l = N;
            }
            else {
                Aini = A+1;
                l=N-1;
            }
            for (size_t i=0; i<R; ++i, Aini += lda+1) {
                if (Q[i] > i){
                    FFLAS::fassign (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
                    for (size_t j=0; j<l-i; ++j)
                        F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
                }
            }
            return R;
        }
    };

    template<>
    class callLUdivine_small<float> {
    public:
        template <class Field>
        inline size_t
        operator()( const Field& F,
                    const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda, size_t*P,
                    size_t *Q, const FFPACK::FFPACK_LU_TAG LuTag)
        {

            if ( !(M && N) ) return 0;
            typedef typename Field::Element elt;
            elt * Aini = A;
            elt * Acurr;
            size_t rowp = 0;
            size_t R = 0;
            size_t k = 0;
            size_t delay =0;
            size_t kmax = FFLAS::Protected::DotProdBoundClassic (F, F.one) -1; // the max number of delayed operations
            while ((rowp<M) && (k<N)){
                size_t colp;

                //Find non zero pivot
                colp = k;
                Acurr = Aini;
                while ((F.isZero(*Acurr)) || (F.isZero (F.reduce (*Acurr))))
                    if (++colp == N){
                        if (rowp==M-1)
                            break;
                        colp=k; ++rowp;
                        Acurr = Aini += lda;
                    }
                    else
                        ++Acurr;

                if ((rowp == M-1)&&(colp == N))
                    break;
                R++;
                P[k] = colp;
                Q[k] = rowp;

                // Permutation of the pivot column
                FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

                //Normalization
                elt invpiv;
                F.init(invpiv);
                F.reduce (*Aini);
                F.inv (invpiv,*Aini);

                for (size_t j=1; j<N-k; ++j)
                    if (!F.isZero(*(Aini+j)))
                        F.reduce (*(Aini+j));
                for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    if (!F.isZero(*(Aini+i)))
                        F.reduce (*(Aini+i));

                if (Diag == FFLAS::FflasUnit) {
                    // for (size_t j=1; j<N-k; ++j)
                    // if (!F.isZero(*(Aini+j)))
                    // F.mulin (*(Aini+j),invpiv);
                    FFLAS::fscalin(F,N-k-1,invpiv,Aini+1,1);
                }
                else {
                    // for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
                    // if (!F.isZero(*(Aini+i)))
                    // F.mulin (*(Aini+i),invpiv);
                    FFLAS::fscalin(F,M-rowp-1,invpiv,Aini+lda,lda);
                }

                if (delay++ >= kmax){ // Reduction has to be done
                    delay = 0;
                    FFLAS::freduce (F, M-rowp-1, N-k-1, Aini+lda+1, lda);
                    // for (size_t i=1; i<M-rowp; ++i)
                    // 	for (size_t j=1; j<N-k; ++j)
                    // 		F.reduce (*(Aini+i*lda+j));
                }
                //Elimination
                for (size_t i=1; i<M-rowp; ++i)
                    for (size_t j=1; j<N-k; ++j)
                        *(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
                //Or equivalently, but without delayed ops :
                //FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

                Aini += lda+1; ++rowp; ++k;
            }

            // Compression the U matrix
            size_t l;
            if (Diag == FFLAS::FflasNonUnit){
                Aini = A;
                l = N;
            }
            else {
                Aini = A+1;
                l=N-1;
            }
            for (size_t i=0; i<R; ++i, Aini += lda+1) {
                if (Q[i] > i){
                    FFLAS::fassign (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
                    for (size_t j=0; j<l-i; ++j)
                        F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
                }
            }
            return R;
        }
    };

    template <class Field>
    inline size_t
    LUdivine (const Field& F,
              const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
              const size_t M, const size_t N,
              typename Field::Element_ptr A, const size_t lda,
              size_t*P, size_t *Q
              , const FFPACK::FFPACK_LU_TAG LuTag // =FFPACK::FfpackSlabRecursive
              , const size_t cutoff // =__FFPACK_LUDIVINE_CUTOFF
             )
    {
        if ( !(M && N) ) return 0;
        typedef typename Field::Element elt;
        size_t MN = std::min(M,N);

        size_t incRow, incCol, rowDim, colDim;
        if (trans == FFLAS::FflasTrans){
            incRow = 1;
            incCol = lda;
            colDim = M;
            rowDim = N;
        }
        else {
            incRow = lda;
            incCol = 1;
            colDim = N;
            rowDim = M;
        }

        if ((rowDim < cutoff) && (colDim < 2*cutoff)) { // the coeff 2 is experimentally determined!
            return LUdivine_small (F, Diag, trans, M, N, A, lda, P, Q, LuTag);
        }
        else { // recursively :
            if (MN == 1){
                size_t ip=0;
                while (F.isZero (*(A+ip*incCol)))
                    if (++ip == colDim)
                        break;
                *Q=0;
                if (ip == colDim){ // current row is zero
                    *P=0;
                    if (colDim == 1){
                        //while (ip<M && !F.isUnit(*(A+ip*lda)))
                        while (ip<rowDim && F.isZero(*(A + ip*incRow))){
                            //	Q[ip]=ip;
                            ip++;
                        }
                        if (ip == rowDim) {
                            return 0;
                        }
                        else {
                            size_t oldip = ip;
                            if ( Diag == FFLAS::FflasNonUnit ){
                                elt invpiv;
                                F.init(invpiv);
                                F.inv(invpiv,*(A+ip*incRow));
                                if (++ip < rowDim)
                                    FFLAS::fscalin(F,rowDim-ip,invpiv,A+ip*incRow,incRow);
                                // elt tmp;
                                //								F.init(tmp);
                                //								F.assign(tmp, *(A+oldip*incRow));
                                //								F.assign( *(A+oldip*incRow), *A);
                                F.assign( *A,*(A+oldip*incRow));
                                F.assign( *(A+oldip*incRow), F.zero);
                            }
                            *Q=oldip;

                            return 1;
                        }
                    }
                    else{
                        *Q=0; return 0;
                    }
                }
                *P=ip;
                if (ip!=0){
                    // swap the pivot
                    typename Field::Element tmp;
                    F.init(tmp);
                    F.assign(tmp,*A);
                    F.assign(*A, *(A + ip*incCol));
                    F.assign(*(A + ip*incCol), tmp);
                }
                elt invpiv;
                F.init(invpiv);
                F.inv(invpiv, *A);
                if ( Diag == FFLAS::FflasUnit && colDim>1){
                    // Normalisation of the row
                    FFLAS::fscalin(F,colDim-1,invpiv,A+incCol,incCol);
                }
                else if ( (colDim==1) &&(Diag==FFLAS::FflasNonUnit) ){
                    if (++ip < rowDim){
                        FFLAS::fscalin(F,rowDim-ip,invpiv,A+ip*incRow,incRow);
                    }
                }
                return 1;
            }
            else { // MN>1
                size_t Nup = rowDim >> 1;
                size_t Ndown =  rowDim - Nup;
                // FFLASFFPACK_check(Ndown < rowDim);
                // Recursive call on NW
                size_t R, R2;
                if (trans == FFLAS::FflasTrans){
                    R = LUdivine (F, Diag, trans, colDim, Nup, A, lda, P, Q,
                                  LuTag, cutoff);

                    typename Field::Element_ptr Ar = A + Nup*incRow;   // SW
                    typename Field::Element_ptr Ac = A + R*incCol;     // NE
                    typename Field::Element_ptr An = Ar+ R*incCol;     // SE

                    if (!R){
                        if (LuTag == FFPACK::FfpackSingular )
                            return 0;
                    }
                    else {
                        FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                                        Ndown, 0,(int) R, Ar, lda, P);
                        // Ar <- L1^-1 Ar
                        FFLAS::ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasLower,
                                      FFLAS::FflasNoTrans, Diag, R, Ndown,
                                      F.one, A, lda, Ar, lda);
                        // An <- An - Ac*Ar
                        if (colDim>R)
                            fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, colDim-R, Ndown, R,
                                   F.mOne, Ac, lda, Ar, lda, F.one, An, lda);
                    }
                    // Recursive call on SE
                    R2 = LUdivine (F, Diag, trans, colDim-R, Ndown, An, lda, P + R, Q + Nup, LuTag, cutoff);
                    for (size_t i = R; i < R + R2; ++i)
                        P[i] += R;
                    if (R2) {
                        // An <- An.P2
                        FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                                        Nup,(int) R, (int)(R+R2), A, lda, P);
                    }
                    else {
                        if (LuTag == FFPACK::FfpackSingular)
                            return 0;
                    }

                }
                else { // trans == FFLAS::FflasNoTrans
                    R = LUdivine (F, Diag, trans, Nup, colDim, A, lda, P, Q, LuTag, cutoff);
                    typename Field::Element_ptr Ar = A + Nup*incRow;   // SW
                    typename Field::Element_ptr Ac = A + R*incCol;     // NE
                    typename Field::Element_ptr An = Ar+ R*incCol;     // SE

                    if (!R){
                        if (LuTag == FFPACK::FfpackSingular )
                            return 0;
                    }
                    else { /*  R>0 */
                        // Ar <- Ar.P^T
                        FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                                        Ndown, 0,(int) R, Ar, lda, P);
                        // Ar <- Ar.U1^-1
                        ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper,
                               FFLAS::FflasNoTrans, Diag, Ndown, R,
                               F.one, A, lda, Ar, lda);
                        // An <- An - Ar*Ac
                        if (colDim>R)
                            fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, colDim-R, R,
                                   F.mOne, Ar, lda, Ac, lda, F.one, An, lda );

                    }
                    // Recursive call on SE
                    R2=LUdivine (F, Diag, trans, Ndown, N-R, An, lda,P+R, Q+Nup, LuTag, cutoff);
                    for (size_t i = R; i < R + R2; ++i)
                        P[i] += R;
                    if (R2)
                        // An <- An.P2
                        FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                                        Nup,(int) R, (int)(R+R2), A, lda, P);
                    else if (LuTag == FFPACK::FfpackSingular)
                        return 0;

                }
                // Non zero row permutations
                for (size_t i = Nup; i < Nup + R2; i++)
                    Q[i] += Nup;
                if (R < Nup){
                    // Permutation of the 0 rows
                    if (Diag == FFLAS::FflasNonUnit){
                        for ( size_t i = Nup, j = R ; i < Nup + R2; ++i, ++j){
                            FFLAS::fassign( F, colDim - j, A + i*incRow + j*incCol, incCol, A + j * (lda + 1), incCol);
                            for (typename Field::Element_ptr Ai = A + i*incRow + j*incCol;
                                 Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
                                F.assign (*Ai, F.zero);
                            ///@todo std::swap ?
                            size_t t = Q[j];
                            Q[j]=Q[i];
                            Q[i] = t;
                        }
                    }
                    else { // Diag == FFLAS::FflasUnit
                        for ( size_t i = Nup, j = R+1 ; i < Nup + R2; ++i, ++j){
                            FFLAS::fassign( F, colDim - j,
                                            A + i*incRow + j*incCol, incCol,
                                            A + (j-1)*incRow + j*incCol, incCol);
                            for (typename Field::Element_ptr Ai = A + i*incRow + j*incCol;
                                 Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
                                F.assign (*Ai, F.zero);
                            size_t t = Q[j-1];
                            Q[j-1]=Q[i];
                            Q[i] = t;
                        }
                    }
                }
                return R + R2;
            }
        }
    }

    namespace Protected {

        //---------------------------------------------------------------------
        // LUdivine_construct: (Specialisation of LUdivine)
        // LUP factorisation of the Krylov base matrix of A^t and v.
        // When all rows have been factorized in A, and rank is full,
        // then new krylov vectors are computed and then triangularized
        // P is the permutation matrix stored in the lapack style
        // nRowX is the number of Krylov vectors already computed,
        // nUsedRowX is the number of Krylov vectors already triangularized
        //---------------------------------------------------------------------

        template <class Field>
        size_t
        LUdivine_construct( const Field& F, const FFLAS::FFLAS_DIAG Diag,
                            const size_t M, const size_t N,
                            typename Field::ConstElement_ptr A, const size_t lda,
                            typename Field::Element_ptr X, const size_t ldx,
                            typename Field::Element_ptr u, const size_t incu, size_t* P,
                            bool computeX
                            , const FFPACK::FFPACK_MINPOLY_TAG MinTag //= FFPACK::FfpackDense
                            , const size_t kg_mc// =0
                            , const size_t kg_mb// =0
                            , const size_t kg_j // =0
                          )
        {

            size_t MN = std::min(M,N);

            if (MN == 1){
                size_t ip=0;
                while (ip<N && F.isZero(*(X+ip))){ip++;}
                if (ip==N){ // current row is zero
                    *P=0;
                    return 0;
                }
                *P=ip;
                if (ip!=0){
                    // swap the pivot
                    F.assign(X[0],X[ip]);
                    F.assign(X[ip],F.zero);
                }
                if ( Diag == FFLAS::FflasUnit ){
                    typename Field::Element invpiv;
                    F.inv(invpiv, *X);

                    // Normalisation of the row
                    // for (size_t k=1; k<N; k++)
                    // F.mulin(*(X+k), invpiv);
                    FFLAS::fscalin(F,N-1,invpiv,X+1,1);
                }
                if (N==1 && M>1 && computeX)// Only appends when A is 1 by 1
                    F.mul(*(X+ldx),*X, *A);

                return 1;
            }
            else{ // MN>1
                size_t Nup = MN>>1;
                size_t Ndown =  M - Nup;

                // Recursive call on NW
                size_t R = LUdivine_construct(F, Diag, Nup, N, A, lda, X, ldx, u, incu,
                                              P, computeX, MinTag, kg_mc, kg_mb, kg_j );
                if (R==Nup){
                    typename Field::Element_ptr Xr = X + Nup*ldx; //  SW
                    typename Field::Element_ptr Xc = X + Nup;     //  NE
                    typename Field::Element_ptr Xn = Xr + Nup;    //  SE
                    typename Field::Element_ptr Xi = Xr;
                    if ( computeX ){
                        if (MinTag == FFPACK::FfpackDense)
                            for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
                                fgemv(F, FFLAS::FflasNoTrans, N, N, F.one,
                                      A, lda, u, incu, F.zero, Xi,1);
                                FFLAS::fassign(F, N,Xi, 1, u,incu);
                            }
                        else // Keller-Gehrig Fast algorithm's matrix
                            for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
                                FFPACK::Protected::fgemv_kgf( F, N, A, lda, u, incu, Xi, 1,
                                                              kg_mc, kg_mb, kg_j );
                                FFLAS::fassign(F, N,Xi, 1, u, incu);
                            }
                    }
                    // Apply the permutation on SW
                    FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
                                    Ndown, 0,(int) R, Xr, ldx, P);
                    // Triangular block inversion of NW and apply to SW
                    // Xr <- Xr.U1^-1
                    ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag,
                           Ndown, R, F.one, X, ldx, Xr, ldx);

                    // Update of SE
                    // Xn <- Xn - Xr*Xc
                    fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, N-Nup, Nup,
                           F.mOne, Xr, ldx, Xc, ldx, F.one, Xn, ldx);

                    // Recursive call on SE

                    size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, A, lda,
                                                   Xn, ldx, u, incu, P + Nup,
                                                   false, MinTag, kg_mc, kg_mb, kg_j);
                    for ( size_t i=R;i<R+R2;++i) P[i] += R;

                    FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
                                    Nup, (int)R, (int)(R+R2), X, ldx, P);

                    return R+=R2;
                }
                else
                    return R;
                // Rank deficient matrices can only be factorized
                // under the condition: the first R rows are linearly independent
                // If not, the lower block is never factorized as soon as the
                // upper block is rank defficient
            }
        }

    } // Protected

} // FFPACK

#endif //__FFLASFFPACK_ffpack_ludivine_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
