/*
 * Copyright (C) 2014 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_ffpack_permutation_INL
#define __FFLASFFPACK_ffpack_permutation_INL


#include <givaro/zring.h>

#include "fflas-ffpack/fflas/fflas_fassign.h"

#define FFLASFFPACK_PERM_BKSIZE 32

namespace FFPACK {
    /** MonotonicApplyP
     * Apply a permutation defined by the first R entries of the vector P (the pivots).
     * The non pivot elements, are located in montonically increasing order.
     */
    template<class Field>
    inline void
    MonotonicApplyP (const Field& F,
                     const FFLAS::FFLAS_SIDE Side,
                     const FFLAS::FFLAS_TRANSPOSE Trans,
                     const size_t M, const size_t ibeg, const size_t iend,
                     typename Field::Element_ptr A, const size_t lda, const size_t * P, const size_t R)
    {
        const size_t B = FFLASFFPACK_PERM_BKSIZE;
        size_t lenP = iend-ibeg;
        size_t * MathP = new size_t[lenP];
        for (size_t i=0; i<lenP; ++i)
            MathP[i] = i;
        LAPACKPerm2MathPerm (MathP, P, lenP);

        std::vector<bool> ispiv(lenP,false);
        size_t pivrowstomove = 0;
        size_t nonpivrowstomove = 0;
        size_t maxpiv = R-1;
        for (size_t i=0; i<R; i++) {
            ispiv[MathP[i]] = true;
            if (MathP[i] != i){
                pivrowstomove++;
                if(maxpiv < MathP[i]) maxpiv = MathP[i];
            }
        }
        if (!pivrowstomove) // Permutation is the identity
            return;

        for (size_t i=R; i<lenP; i++)
            if (MathP[i] != i)
                nonpivrowstomove++;
        size_t NB = M/B;
        size_t last = M%B;
        size_t incA, llda;
        if (Side == FFLAS::FflasLeft)  {incA = 1; llda = lda;}
        else {incA = lda; llda = 1;}
        size_t inc = B*incA;

        if (((Side == FFLAS::FflasLeft) && (Trans == FFLAS::FflasNoTrans)) ||
            ((Side == FFLAS::FflasRight) && (Trans == FFLAS::FflasTrans))){
            // Compressing
#ifdef MONOTONIC_CYCLES
            for (size_t i = 0; i<NB; i++)
                MonotonicCompressCycles (F, Side, B, A+i*inc, llda, incA, MathP, lenP);
            MonotonicCompressCycles (F, Side, last, A+NB*inc, llda, incA, MathP, lenP);
#elif defined MONOTONIC_MOREPIVOTS
            for (size_t i = 0; i<NB; i++)
                MonotonicCompressMorePivots (F, Side, B, A+i*inc, llda, incA, MathP, R, nonpivrowstomove, lenP);
            MonotonicCompressMorePivots (F, Side, last, A+NB*inc, llda, incA, MathP, R, nonpivrowstomove, lenP);
#else
            for (size_t i = 0; i<NB; i++)
                MonotonicCompress (F, Side, B, A+i*inc, llda, incA, MathP, R, maxpiv, pivrowstomove, ispiv);
            MonotonicCompress (F, Side, last, A+NB*inc, llda, incA, MathP, R, maxpiv, pivrowstomove, ispiv);
#endif
        } else {
            // Expanding
            for (size_t i = 0; i<NB; i++)
                MonotonicExpand (F, Side, B, A+i*inc, llda, incA, MathP, R, maxpiv, pivrowstomove, ispiv);
            MonotonicExpand (F, Side, last, A+NB*inc, llda, incA, MathP, R, maxpiv, pivrowstomove, ispiv);
        }
        delete[] MathP;
    }

    template<class Field>
    inline void
    MonotonicCompress (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                       typename Field::Element_ptr A, const size_t lda, const size_t incA,
                       const size_t * MathP, const size_t R, const size_t maxpiv,
                       const size_t rowstomove, const std::vector<bool> &ispiv)
    {
        // Storing pivot rows in temp
        typename Field::Element_ptr temp= FFLAS::fflas_new (F, rowstomove, M);
        size_t ldtemp=M;
        for (size_t i=0,j=0; i<R; i++){
            if (MathP[i] != i){
                FFLAS::fassign (F, M, A+MathP[i]*lda, incA, temp+j*ldtemp, 1);
                j++;
            }
        }
        // Moving non pivot rows to the R+1 .. iend positions
        int dest = maxpiv;
        int src = dest - 1;
        while (dest >= (int)R){
            if ((src >= 0) && ispiv[src]){ // src points to a pivot row: skip it
                src--;
                continue;
            }
            FFLAS::fassign(F, M, A+src*lda, incA, A+dest*lda, incA);
            src--; dest--;
        }
        // Moving the pivots to their position in the first R rows
        for (size_t i=0, j=0; i<R; i++)
            if (MathP[i] != i){
                FFLAS::fassign (F, M, temp + j*ldtemp, 1, A + i*lda, incA);
                j++;
            }
        FFLAS::fflas_delete(temp);
    }

    template<class Field>
    inline void
    MonotonicCompressMorePivots (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                                 typename Field::Element_ptr A, const size_t lda, const size_t incA,
                                 const size_t * MathP, const size_t R, const size_t rowstomove, const size_t lenP)
    {
        std::vector<bool> done(lenP,false);
        typename Field::Element_ptr temp= FFLAS::fflas_new (F, rowstomove, M);
        size_t ldtemp=M;
        // Move every non pivot row to temp
#ifdef VERBOSE
        std::cerr<<"R = "<<R<<std::endl;
        write_perm(std::cerr<<"MathP = ",MathP,lenP);
#endif
        for (size_t i=R,j=0; i<lenP; i++){
            if (MathP[i] != i){
#ifdef VERBOSE
                std::cerr<<"A["<<MathP[i]<<"] -> temp["<<j<<"]"<<std::endl;
#endif
                FFLAS::fassign (F, M, A+MathP[i]*lda, incA, temp+j*ldtemp, 1);
                done[MathP[i]]=true;
                j++;
            }
        }
        // Move the pivot rows of every cycle containing a non pivot row (avoiding to use a temp)
        for (size_t i=R; i<lenP; i++){
            size_t j=MathP[i];
            while ((MathP[j] != j) && (!done[MathP[j]])){
                // A[P[j]] -> A[j]
#ifdef VERBOSE
                std::cerr<<"Moving pivots 1 A["<<MathP[j]<<"] -> A["<<j<<"]"<<std::endl;
#endif
                FFLAS::fassign (F, M, A+MathP[j]*lda, incA, A+j*lda, incA);
                done[MathP[j]] = true;
                j = MathP[j];
            }
        }

        // Moving the remaining cycles using one vector temp
        typename Field::Element_ptr tmprow = FFLAS::fflas_new(F,M);
        for (size_t i=0; i<R; i++){
            if ((MathP[i]!=i)&&(!done[MathP[i]])){ // entering a cycle
                size_t j=i;
#ifdef VERBOSE
                std::cerr<<"Moving pivots 2 A["<<j<<"] -> tmprow"<<std::endl;
#endif
                FFLAS::fassign (F, M, A+j*lda, incA, tmprow, 1);
                done[j] = true;
                do{
                    // A[P[j]] -> A[j]
#ifdef VERBOSE
                    std::cerr<<"Moving pivots 2 A["<<MathP[j]<<"] -> A["<<j<<"]"<<std::endl;
#endif
                    FFLAS::fassign (F, M, A+MathP[j]*lda, incA, A+j*lda, incA);
                    done[MathP[j]] = true;
                    j = MathP[j];
                } while (!done[MathP[j]]);
                FFLAS::fassign (F, M, tmprow, 1, A+j*lda, incA);
#ifdef VERBOSE
                std::cerr<<"Moving pivots 2 tmprow -> A["<<j<<"]"<<std::endl;
#endif
            }
        }
        // Move the non pivot rows to the last lenP-R positions
        for (size_t i=R, j=0; i<lenP; i++)
            if (MathP[i] != i){
#ifdef VERBOSE
                std::cerr<<"temp["<<j<<"] -> A["<<i<<"] "<<std::endl;
#endif
                FFLAS::fassign (F, M, temp + j*ldtemp, 1, A + i*lda, incA);
                j++;
            }

        FFLAS::fflas_delete(tmprow);
        FFLAS::fflas_delete(temp);
    }

    template<class Field>
    inline void
    MonotonicCompressCycles (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                             typename Field::Element_ptr A, const size_t lda, const size_t incA,
                             const size_t * MathP, const size_t lenP)
    {
        std::vector<bool> done(lenP,false);
        // Move every non pivot row to temp
#ifdef VERBOSE
        write_perm(std::cerr<<"MathP = ",MathP,lenP);
#endif
        // Moving the remaining cycles using one vector temp
        typename Field::Element_ptr tmprow = FFLAS::fflas_new(F,FFLASFFPACK_PERM_BKSIZE);
        for (size_t i=0; i<lenP; i++){
            if ((MathP[i]!=i)&&(!done[MathP[i]])){ // entering a cycle
                size_t j=i;
#ifdef VERBOSE
                std::cerr<<"Moving pivots A["<<j<<"] -> tmprow"<<std::endl;
#endif
                FFLAS::fassign (F, M, A+j*lda, incA, tmprow, 1);
                done[j] = true;
                do{
                    // A[P[j]] -> A[j]
#ifdef VERBOSE
                    std::cerr<<"Moving pivots A["<<MathP[j]<<"] -> A["<<j<<"]"<<std::endl;
#endif
                    FFLAS::fassign (F, M, A+MathP[j]*lda, incA, A+j*lda, incA);
                    done[MathP[j]] = true;
                    j = MathP[j];
                } while (!done[MathP[j]]);
                FFLAS::fassign (F, M, tmprow, 1, A+j*lda, incA);
#ifdef VERBOSE
                std::cerr<<"Moving pivots tmprow -> A["<<j<<"]"<<std::endl;
#endif
            }
        }
        FFLAS::fflas_delete(tmprow);
    }
    template<class Field>
    void
    MonotonicExpand (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                     typename Field::Element_ptr A, const size_t lda, const size_t incA,
                     const size_t * MathP, const size_t R, const size_t maxpiv,
                     const size_t rowstomove, const std::vector<bool> &ispiv)
    {
        // Storing pivot rows in temp
        typename Field::Element_ptr temp= FFLAS::fflas_new (F, rowstomove, M);
        size_t ldtemp=M;
        for (size_t i=0,j=0; i<R; i++){
            if (MathP[i] != i){
                FFLAS::fassign (F, M, A+i*lda, incA, temp+j*ldtemp, 1);
                j++;
            }
        }
        // Moving the non pivot rows
        size_t dest = 0;
        size_t src = R;
        while (src <= maxpiv){
            if (ispiv[dest]){ // src points to a pivot row: skip it
                dest++;
                continue;
            }
            FFLAS::fassign(F, M, A+src*lda, incA, A+dest*lda, incA);
            src++; dest++;
        }
        // Moving the pivots to their final position
        for (size_t i=0, j=0; i<R; i++)
            if (MathP[i] != i){
                FFLAS::fassign (F, M, temp + j*ldtemp, 1, A + MathP[i]*lda, incA);
                j++;
            }
        FFLAS::fflas_delete(temp);
    }

    template<class Field>
    inline void
    applyP_block (const Field& F,
                  const FFLAS::FFLAS_SIDE Side,
                  const FFLAS::FFLAS_TRANSPOSE Trans,
                  const size_t M, const size_t ibeg, const size_t iend,
                  typename Field::Element_ptr A, const size_t lda, const size_t * P)
    {
        if ( Side == FFLAS::FflasRight ) {
            if ( Trans == FFLAS::FflasTrans ){
                for ( size_t i=(size_t)ibeg; i<(size_t) iend; ++i)
                    if ( P[i]!= i )
                        FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda);
            } else { // Trans == FFLAS::FflasNoTrans
                for (size_t i=iend; i-->ibeg; )
                    if ( P[i]!=(size_t)i )
                        FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda);
            }
        } else { // Side == FFLAS::FflasLeft
            if ( Trans == FFLAS::FflasNoTrans ) {
                for (size_t i=(size_t)ibeg; i<(size_t)iend; ++i)
                    if ( P[i]!= (size_t) i )
                        FFLAS::fswap( F, M, A + P[i]*lda, 1, A + i*lda, 1);
            } else { // Trans == FFLAS::FflasTrans
                for (size_t i=iend; i-->ibeg; )
                    if ( P[i]!= (size_t) i )
                        FFLAS::fswap( F, M, A + P[i]*lda, 1, A + i*lda, 1);
            }
        }
    }

    template<class Field>
    inline void applyP( const Field& F,
                        const FFLAS::FFLAS_SIDE Side,
                        const FFLAS::FFLAS_TRANSPOSE Trans,
                        const size_t M, const size_t ibeg, const size_t iend,
                        typename Field::Element_ptr A, const size_t lda, const size_t * P,
                        const FFLAS::ParSeqHelper::Sequential seq)
    {

        const size_t bk = FFLASFFPACK_PERM_BKSIZE;
        const size_t NB = M/bk;
        const size_t last = M%bk;
        const size_t incA = (Side == FFLAS::FflasLeft)? 1:lda;
        const size_t inc = bk*incA;

        for (size_t i = 0; i<NB; i++)
            applyP_block (F, Side, Trans, bk, ibeg, iend, A+i*inc, lda, P);
        applyP_block (F, Side, Trans, last, ibeg, iend, A+NB*inc, lda, P);
    }

    template<class Field>
    inline void doApplyS (const Field& F,
                          typename Field::Element_ptr A, const size_t lda, typename Field::Element_ptr tmp,
                          const size_t width, const size_t M2,
                          const size_t R1, const size_t R2,
                          const size_t R3, const size_t R4)
    {
        FFLAS::fassign(F, M2-R1-R2, width,  A + (R1+R2)*lda, lda, tmp, width);
        FFLAS::fassign(F, R3+R4, width,  A + M2*lda, lda, A + (R1+R2)*lda, lda);
        FFLAS::fassign(F, M2-R1-R2, width, tmp, width, A + (R1+R2+R3+R4)*lda, lda);
    }
    template <class Field>
    inline void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
                              const size_t width, const size_t M2,
                              const size_t R1, const size_t R2,
                              const size_t R3, const size_t R4){
        MatrixApplyS (F, A, lda, width, M2, R1, R2, R3, R4, FFLAS::ParSeqHelper::Sequential());
    }
    template <class Field>
    inline void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
                              const size_t width, const size_t M2,
                              const size_t R1, const size_t R2,
                              const size_t R3, const size_t R4,
                              const FFLAS::ParSeqHelper::Sequential seq)
    {
        typename Field::Element_ptr tmp = FFLAS::fflas_new (F, M2-R1-R2, width);
        doApplyS (F, A, lda, tmp, width, M2, R1, R2, R3, R4);
        FFLAS::fflas_delete (tmp);
    }
    template <class Field, class Cut, class Param>
    inline void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
                              const size_t width, const size_t M2,
                              const size_t R1, const size_t R2,
                              const size_t R3, const size_t R4,
                              const FFLAS::ParSeqHelper::Parallel<Cut, Param> par)
    {
        SYNCH_GROUP(
            FORBLOCK1D(iter,width, par,
                       TASK(MODE(CONSTREFERENCE(F,A) READ(A[BLOCKSIZE*t])),
                            MatrixApplyS (F, A+iter.begin(), lda, iter.end()-iter.begin(), M2, R1, R2, R3, R4););
                       );
                   );
    }

    template <class T>
    inline void PermApplyS (T* A, const size_t lda,
                            const size_t width, const size_t M2,
                            const size_t R1, const size_t R2,
                            const size_t R3, const size_t R4)
    {
        Givaro::ZRing<T> D;
        T* tmp = FFLAS::fflas_new<T>((M2-R1-R2)*width);
        doApplyS (D, A, lda, tmp, width, M2, R1, R2, R3, R4);
        FFLAS::fflas_delete( tmp);
    }



    template <class Field>
    inline void doApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda, typename Field::Element_ptr tmp,
                          const size_t width, const size_t N2,
                          const size_t R1, const size_t R2,
                          const size_t R3, const size_t R4)
    {
        for (size_t k = 0; k < width; ++k){
            FFLAS::fassign(F, N2-R1, A+R1+k*lda, 1, tmp + k*(N2-R1), 1);
            FFLAS::fassign(F, R2, A+N2+k*lda, 1, A + R1 + k*lda, 1);
            FFLAS::fassign(F, R3, tmp + k*(N2-R1), 1, A+R1+R2+k*lda, 1);
            FFLAS::fassign(F, R4, A + N2 + R2 + k*lda, 1, A + R1+R2+R3 + k*lda, 1);
            FFLAS::fassign(F, N2-R1-R3, tmp + R3 + k*(N2-R1), 1, A+R1+R2+R3+R4+k*lda, 1);
        }
    }

    template <class Field>
    inline void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
                              const size_t width, const size_t N2,
                              const size_t R1, const size_t R2,
                              const size_t R3, const size_t R4){
        MatrixApplyT (F, A, lda, width, N2, R1, R2, R3, R4, FFLAS::ParSeqHelper::Sequential());
    }

    template <class Field>
    inline void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
                              const size_t width, const size_t N2,
                              const size_t R1, const size_t R2,
                              const size_t R3, const size_t R4,
                              const FFLAS::ParSeqHelper::Sequential seq)
    {
        typename Field::Element_ptr tmp = FFLAS::fflas_new (F, N2-R1, width);
        doApplyT (F, A, lda, tmp, width, N2, R1, R2, R3, R4);
        FFLAS::fflas_delete (tmp);
    }

    template <class Field, class Cut, class Param>
    inline void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
                       const size_t width, const size_t N2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4,
                       const FFLAS::ParSeqHelper::Parallel<Cut, Param> par)
    {
        SYNCH_GROUP(
            FORBLOCK1D(iter, width, par,
                       TASK(MODE(CONSTREFERENCE(F, A) READWRITE(A[BLOCKSIZE*t*lda])),
                            MatrixApplyT(F,A+iter.begin()*lda, lda, iter.end()-iter.begin(), N2, R1, R2, R3, R4) );
                       );
                    );
    }

    template <class T>
    inline void PermApplyT (T* A, const size_t lda,
                            const size_t width, const size_t N2,
                            const size_t R1, const size_t R2,
                            const size_t R3, const size_t R4)
    {
        Givaro::ZRing<T> D;
        T* tmp = FFLAS::fflas_new<T >((N2-R1)*width);
        doApplyT (D, A, lda, tmp, width, N2, R1, R2, R3, R4);
        FFLAS::fflas_delete( tmp);
    }

    /**
     * Conversion of a permutation from LAPACK format to Math format
     */
    inline void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
                                     const size_t N)
    {
        for (size_t i=0; i<N; i++)
            MathP[i] = i;
        for (size_t i=0; i<N; i++){
            if (LapackP[i] != i){
                std::swap(MathP[i],MathP[LapackP[i]]);
            }
        }
    }

    /**
     * Conversion of a permutation from Maths format to LAPACK format
     */
    inline void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
                                     const size_t N)
    {
        size_t * T = FFLAS::fflas_new<size_t>(N);
        size_t * Tinv = FFLAS::fflas_new<size_t>(N);
        for (size_t i=0; i<N; i++){
            T[i] =i;
            Tinv[i] = i;
        }
        for (size_t i=0; i<N; i++){
            size_t j = Tinv [MathP [i]];
            LapackP [i] = j;
            size_t tmp = T[j];
            T[j] = T[i];
            Tinv[T[i]] = j;
            T[i] = tmp;
            Tinv[tmp] = i;
        }
        FFLAS::fflas_delete( T);
        FFLAS::fflas_delete( Tinv);
    }

    /**
     * @brief Computes P1 x Diag (I_R, P2) where P1 is a LAPACK and P2 a LAPACK permutation
     * and store the result in P1 as a LAPACK permutation
     * @param [inout] P1 a LAPACK permutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */
    inline void composePermutationsLLL (size_t * P1,
                                        const size_t * P2,
                                        const size_t R, const size_t N){
        size_t * MathP = FFLAS::fflas_new<size_t >(N);
        composePermutationsLLM(MathP, P1, P2, R, N);
        MathPerm2LAPACKPerm (P1, MathP, N);
        FFLAS::fflas_delete (MathP);
    }
    /**
     * @brief Computes P1 x Diag (I_R, P2) where P1 is a LAPACK and P2 a LAPACK permutation
     * and store the result in MathP as a MathPermutation format.
     * @param [out]MathP a MathPermutation of size N
     * @param P1 a LAPACK permutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */
    inline void composePermutationsLLM (size_t * MathP,
                                        const size_t * P1,
                                        const size_t * P2,
                                        const size_t R, const size_t N)
    {
        for (size_t i=0; i<N; ++i)
            MathP[i] = i;
        LAPACKPerm2MathPerm (MathP, P1, N);
        composePermutationsMLM (MathP, P2, R, N);
    }

    /**
     * @brief Computes MathP1 x Diag (I_R, P2) where MathP1 is a MathPermutation and P2 a LAPACK permutation
     * and store the result in MathP1 as a MathPermutation format.
     * @param [inout] MathP1 a MathPermutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */
    inline void composePermutationsMLM (size_t * MathP1,
                                        const size_t * P2,
                                        const size_t R, const size_t N)
    {
        for (size_t i=R; i<N; i++){
            if (P2[i-R] != i-R){
                size_t tmp = MathP1[i];
                MathP1[i] = MathP1[P2[i-R]+R];
                MathP1[P2[i-R]+R] = tmp;
            }
        }
    }

    inline void
    cyclic_shift_mathPerm (size_t * P,  const size_t s)
    {
        size_t tmp;
        tmp = P[s-1];
        //memmove(P+1, P, (s)*sizeof(size_t));
        size_t * Pi = P;
        std::copy(Pi, Pi+s-1, Pi+1);

        *(P)=tmp;
    }
    // @BUG highly not portable to other fields than modular<basis type>
    // Need a rewrite in order to support RNSModP field
    template<class Field>
    inline void cyclic_shift_row_col(const Field & F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda)
    {
        typedef typename Field::Element Element;
        typedef typename Field::Element_ptr Element_ptr;
#ifdef MEMCOPY
        //		std::cerr << "BEF m: " << m << ", n: " << n << std::endl;

        if (m > 1) {
            const size_t mun(m-1);
            if (n > 1) {
                //     std::cerr << "m: " << m << ", n: " << n << std::endl;
                const size_t nun(n-1);
                const size_t blo(sizeof(Element));
                // const size_t bmu(blo*mun);
                const size_t bnu(blo*nun);
                Element_ptr b = FFLAS::fflas_new(F,mun);
                for(size_t i=0; i<mun; ++i) b[i] = A[i*lda+nun];
                Element_ptr dc = FFLAS::fflas_new (F,n);
                memcpy(dc+1,A+mun*lda,bnu);
                *dc = *(A+mun*lda+nun); // this is d
                // dc = [ d c ]

                for(size_t i=mun; i>0; --i)
                    memcpy(A+1+i*lda, A+(i-1)*lda, bnu);

                memcpy(A, dc, bnu+blo);
                for(size_t i=0; i<mun; ++i) A[(i+1)*lda] = b[i];
                delete [] dc;
                delete [] b;

            } else if (n != 0) {
                Base_t d = A[mun*lda];
                for(size_t i=mun; i>0; --i) A[i*lda]=A[(i-1)*lda];
                *A=d;
            }
        } else {
            if ((m!=0) && (n > 1)) {
                const size_t nun(n-1);
                const size_t blo(sizeof(Element));
                const size_t bnu(blo*nun);
                Element d = A[nun];
                //  std::cerr << "d: " << d << std::endl;
                Element_ptr tmp = FFLAS::fflas_new(F,nun);
                memcpy(tmp,A,bnu);
                memcpy(A+1,tmp,bnu);
                //				std::copy(A,A+nun,A+1);
                *A=d;
                delete [] tmp;
            }
        }
        //		std::cerr << "AFT m: " << m << ", n: " << n << std::endl;

#else

        //	std::cerr << "BEF m: " << m << ", n: " << n << std::endl;
        if (m > 1) {
            const size_t mun(m-1);
            if (n > 1) {
                const size_t nun(n-1);

                Element_ptr b = FFLAS::fflas_new (F,mun);
                Element_ptr Ainun = A+nun;
                for(size_t i=0; i<mun; ++i, Ainun+=lda) b[i] = *Ainun;

                // dc = [ d c ]
                Element_ptr dc = FFLAS::fflas_new (F,n);
                FFLAS::fassign(F,nun,Ainun-nun,1, dc+1,1);
                //std::copy(Ainun-nun, Ainun, dc+1);

                // this is d
                *dc = *Ainun;

                Element_ptr Ai = A+(mun-1)*lda;
                for(size_t i=mun; i>0; --i, Ai-=lda)
                    FFLAS::fassign(F, nun, Ai,1,Ai+1+lda,1);
                //				std::copy(Ai, Ai+nun, Ai+1+lda);

                FFLAS::fassign(F, n, dc, 1, A, 1);
                //std::copy(dc, dc+n, A);

                Element_ptr Aipo = A+lda;
                for(size_t i=0; i<mun; ++i, Aipo+=lda) *Aipo = b[i];

                FFLAS::fflas_delete(dc);
                FFLAS::fflas_delete(b);
            } else if (n != 0) {
                Element_ptr Ai=A+mun*lda;
                Element_ptr d = *Ai;
                for(; Ai != A; Ai-=lda) *Ai= *(Ai-lda);
                *A=d;
            }
        } else {
            if ((m!=0) && (n > 1)) {
                const size_t nun(n-1);
                Element d = A[nun];
                FFLAS::fassign(F,nun,A,1,A+1,1);
                //std::copy(A,A+nun,A+1);
                *A=d;
            }
        }

#endif
    }

    template<class Field>
    inline void cyclic_shift_row(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda)
    {

#ifdef MEMCOPY
        if (m > 1) {
            const size_t mun(m-1);

            typename Field::Element_ptr b = FFLAS::fflas_new (F,n);
            typename Field::Element_ptr Ai = A+mun*lda;
            //@BUG not safe with RNSModp field
            memcpy (b,Ai,n*sizeof(typename Field::Element));

            for(typename Field::Element_ptr Ac = A+mun*lda; Ac!=A;Ac-=lda)
                memcpy (Ac, Ac-lda, n*sizeof(typename Field::Element));

            memcpy ( A, b, n*sizeof(typename Field::Element));
            FFLAS::fflas_delete (b);
        }

#else
        if (m > 1) {
            const size_t mun(m-1);

            typename Field::Element_ptr b = FFLAS::fflas_new (F, n);
            typename Field::Element_ptr Ai = A+mun*lda;
            for(size_t i=0; i<n; ++i, Ai+=1) b[i] = *Ai;

            for(typename Field::Element_ptr Ac = A+mun*lda; Ac!=A;Ac-=lda)
                FFLAS::fassign(F,n, Ac-lda, 1, Ac, 1);
            //std::copy(Ac-lda,Ac-lda+n, Ac);

            typename Field::Element_ptr Aii = A;
            for(size_t i=0; i<n; ++i, Aii+=1) *Aii = b[i];

            FFLAS::fflas_delete (b);
        }

#endif
    }

    template<typename T>
    inline void cyclic_shift_row(const RNSIntegerMod<T>& F, typename T::Element_ptr A, size_t m, size_t n, size_t lda)
    {
        if (m > 1) {
            const size_t mun(m-1);

            typename T::Element_ptr b = FFLAS::fflas_new (F, n, 1);
            typename T::Element_ptr Ai = A+mun*lda;
            for(size_t i=0; i<n; ++i, Ai+=1) F.assign(b[i] , *Ai);

            for(typename T::Element_ptr Ac = A+mun*lda; Ac!=A;Ac-=lda)
                FFLAS::fassign(F, n, Ac-lda, 1, Ac, 1);

            typename T::Element_ptr Aii = A;
            for(size_t i=0; i<n; ++i, Aii+=1) F.assign(*Aii, b[i]);

            FFLAS::fflas_delete (b);
        }
    }

    template<class Field>
    inline void cyclic_shift_col(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda)
    {
        if (n > 1) {
            const size_t nun(n-1);
            for(typename Field::Element_ptr Ai=A; Ai!= A+m*lda; Ai+=lda)
            {
                typename Field::Element tmp;
                F.init(tmp);
                F.assign(tmp, Ai[nun]);
                //@BUG: not safe with RNSModP field
                std::copy_backward(Ai, Ai+nun, Ai+n);
                *Ai=tmp;
            }
        }
    }

    template<typename T>
    inline void cyclic_shift_col(const RNSIntegerMod<T>& F, typename T::Element_ptr A, size_t m, size_t n, size_t lda)
    {
        if (n > 1) {
            const size_t nun(n-1);
            for(typename T::Element_ptr Ai=A; Ai!= A+m*lda; Ai+=lda)
            {
                typename T::Element tmp; F.init(tmp);
                F.assign(tmp, Ai[nun]);
                //std::copy_backward(Ai, Ai+nun, Ai+n);
                typename T::Element_ptr Xi = Ai+nun;
                typename T::ConstElement_ptr Yi=Ai+nun-1;
                for (size_t i =0;i<nun;++i, --Xi, --Yi)
                    F.assign(*Xi,*Yi);
                F.assign(*Ai,tmp);
            }
        }
    }

    template<class Field>
    inline void applyP( const Field& F,
                        const FFLAS::FFLAS_SIDE Side,
                        const FFLAS::FFLAS_TRANSPOSE Trans,
                        const size_t m, const size_t ibeg, const size_t iend,
                        typename Field::Element_ptr A, const size_t lda, const size_t * P)
    {
        applyP(F, Side, Trans, m, ibeg, iend, A, lda, P, FFLAS::ParSeqHelper::Sequential());
    }
    //#if defined(__FFLASFFPACK_USE_OPENMP) and defined(_OPENMP)
    template<class Field, class Cut, class Param>
    inline void applyP( const Field& F,
                        const FFLAS::FFLAS_SIDE Side,
                        const FFLAS::FFLAS_TRANSPOSE Trans,
                        const size_t m, const size_t ibeg, const size_t iend,
                        typename Field::Element_ptr A, const size_t lda, const size_t * P,
                        const FFLAS::ParSeqHelper::Parallel<Cut, Param> PSH)
    {
        size_t incBK = (Side == FFLAS::FflasRight)?lda:1;
        SYNCH_GROUP(
            FORBLOCK1D(iter, m, PSH,
                       TASK(MODE(CONSTREFERENCE(F, A,P) READWRITE(A[iter.begin()*incBK])),
                            applyP(F, Side, Trans, iter.end()-iter.begin(), ibeg, iend, A+iter.begin()*incBK, lda, P));
                       );
                    );
    }

} // FFPACK

#endif // __FFLASFFPACK_ffpack_permutation_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
