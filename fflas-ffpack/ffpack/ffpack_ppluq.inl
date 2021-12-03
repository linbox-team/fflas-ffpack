/* ffpack/ffpack_ppluq.inl
 * Copyright (C) 2014 Ziad Sultan
 *
 * Written by Ziad.Sultan@imag.fr
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
 * License along with this library; if not, WRITE to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_ffpack_ppluq_INL
#define __FFLASFFPACK_ffpack_ppluq_INL


//#ifdef __FFLASFFPACK_USE_OPENMP

#define __FFLAS__TRSM_READONLY

#define PBASECASE_K 256


namespace FFPACK {

    template<class Field>
    void  threads_fgemm(const size_t m, const size_t n, const size_t r, int nbthreads, size_t * W1, size_t * W2, size_t * W3, size_t gamma)
    {
        size_t H1, H2, H3;
        size_t M2 = m>>1;
        size_t N2 = n>>1;

        H1 = ((m-N2)*r*(N2-r))<<1;
        H2 = ((M2-r)*r*(n-N2))<<1;
        H3 = ((m-M2)*r*(n-N2))<<1;

        // if we take into account 2 concurrent pluq calls....
        size_t h;
        size_t z1= h*((m-M2)*(N2-r)*(N2-r)-(N2-r)*(N2-r)*(N2-r)/3);
        size_t z2= h*((n-N2)*(M2-r)*(M2-r)-(M2-r)*(M2-r)*(M2-r)/3);

        H1+= z1;
        H2+= z2;

        // compute number of threads for each fgemm call
        *W1=std::max(H1*nbthreads/(H1+H2+H3),(size_t)1);
        *W2=std::max(H2*nbthreads/(H1+H2+H3),(size_t)1);
        *W3=std::max(nbthreads-*W1-*W2,(size_t)1);

        // add gamma factor to change number of threads for pluq calls
        W1-= gamma*z1/(z1+z2);
        W2-= gamma*(1-z1/(z1+z2));
        W3+= gamma;

    }

    template<class Field>
    void threads_ftrsm(const size_t m, const size_t n, int nbthreads, size_t * t1, size_t * t2)
    {
        *t1 = nbthreads*m/(m+n);
        *t2 = nbthreads-(int)*t1;
    }


    // TODO: replace pPLUQ and "int nt", by PLUQ and a Parallel Helper ...
    template<class Field>
    inline size_t
    PLUQ (const Field& Fi, const FFLAS::FFLAS_DIAG Diag, const size_t M, const size_t N,
          typename Field::Element_ptr A, const size_t lda, size_t* P, size_t* Q,
          const FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                              FFLAS::StrategyParameter::Threads>& PSHelper)
    {

        for (size_t i=0; i<M; ++i) P[i] = i;
        for (size_t i=0; i<N; ++i) Q[i] = i;
        if (std::min(M,N) == 0) return 0;
        if (std::max (M,N) == 1) return (Fi.isZero(*A))? 0 : 1;
        if (M == 1){
            size_t piv = 0;
            while ((piv < N) && Fi.isZero (A[piv])) piv++;
            if (piv == N)
                return 0;
            if (piv){
                Q[0] = piv;
                Fi.assign (*A, A[piv]);
                Fi.assign (A[piv], Fi.zero);
            }
            if (Diag== FFLAS::FflasUnit){
                typename Field::Element invpivot;
                Fi.inv(invpivot, *A);
                // for (size_t i=piv+1; i<N; ++i)
                // Fi.mulin (A[i], invpivot);
                FFLAS::fscalin(Fi,N-piv-1,invpivot,A+piv+1,1);
            }
            return 1;
        }
        if (N == 1){
            size_t piv = 0;
            while ((piv < M) && Fi.isZero (A[piv*lda])) piv++;
            if (piv == M)
                return 0;
            if (piv){
                P[0] = piv;
                Fi.assign (*A, *(A+piv*lda));
                Fi.assign (*(A+piv*lda), Fi.zero);
            }
            if (Diag== FFLAS::FflasNonUnit){
                typename Field::Element invpivot;
                Fi.inv(invpivot, *A);
                // for (size_t i=piv+1; i<M; ++i)
                // Fi.mulin (*(A+i*lda), invpivot);
                FFLAS::fscalin(Fi,M-piv-1,invpivot,A+(piv+1)*lda,lda);
            }
            return 1;
        }

#ifdef __FFLASFFPACK_PLUQ_THRESHOLD
        if (std::min(M,N) < __FFLASFFPACK_PLUQ_THRESHOLD)
            return PLUQ_basecaseCrout (Fi, Diag, M, N, A, lda, P, Q);
#endif
        FFLAS::FFLAS_DIAG OppDiag = (Diag == FFLAS::FflasUnit)? FFLAS::FflasNonUnit : FFLAS::FflasUnit;

        size_t M2 = M >> 1;
        size_t N2 = N >> 1;
        size_t * P1 = FFLAS::fflas_new<size_t> (M2);
        size_t * Q1 = FFLAS::fflas_new<size_t> (N2);
        size_t* MathP = 0;
        size_t* MathQ = 0;
        size_t* P2,*P3,*Q2,*Q3,*P4,*Q4;
        size_t R1,R2,R3,R4;

        // A1 = P1 [ L1 ] [ U1 V1 ] Q1
        //        [ M1 ]
        R1 = PLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1, PSHelper);

        typename Field::Element * A2 = A + N2;
        typename Field::Element * A3 = A + M2*lda;
        typename Field::Element * A4 = A3 + N2;
        typename Field::Element * F = A2 + R1*lda;
        typename Field::Element * G = A3 + R1;

        // const FFLAS::CuttingStrategy meth = FFLAS::RECURSIVE;
        // const FFLAS::StrategyParameter strat = FFLAS::TWO_D_ADAPT;

        int nt = PSHelper.numthreads();
        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                               FFLAS::StrategyParameter::TwoDAdaptive> MMParH (std::max(nt/2,1));
        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                               FFLAS::StrategyParameter::TwoDAdaptive> MMFullParH (std::max(nt,1));
        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,
                                               FFLAS::StrategyParameter::Threads> TRSMParH (std::max(nt/2,1));
        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                               FFLAS::StrategyParameter::Threads> PLUQParH1 (std::max(nt/2,1));
        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                               FFLAS::StrategyParameter::Threads> PLUQParH2 (std::max(nt-nt/2,1));

        typename FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,
                                               FFLAS::StrategyParameter::Threads> PermParH (std::max(nt/2,1));
        SYNCH_GROUP(

                    // [ B1 ] <- P1^T A2
                    // [ B2 ]
                    TASK(MODE(READ(P1) CONSTREFERENCE(Fi, P1, A2) READWRITE(A2[0])),
                         { applyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M2, A2, lda, P1, PermParH); }
                        );
                    // [ C1 C2 ] <- A3 Q1^T
                    TASK(MODE(READ(Q1) CONSTREFERENCE(Fi, Q1, A3) READWRITE(A3[0])),
                         applyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N2, A3, lda, Q1,PermParH););
                    // [ C1 C2 ] <- A3 Q1^T

                    CHECK_DEPENDENCIES;
                    // D <- L1^-1 B1
                    TASK(MODE(READ(A[0], R1) CONSTREFERENCE(Fi, TRSMParH, A2) READWRITE(A2[0])),
                         ftrsm( Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R1, N-N2, Fi.one, A, lda, A2, lda, TRSMParH));

                    // E <- C1 U1^-1
                    TASK(MODE(READ(R1, A[0]) CONSTREFERENCE(A3, Fi, M2, R1, TRSMParH) READWRITE(A3[0])),
                         ftrsm(Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R1, Fi.one, A, lda, A3, lda,  TRSMParH));

                    CHECK_DEPENDENCIES;

                    // F <- B2 - M1 D
                    TASK(MODE(READ(A2[0], A[R1*lda], MMParH) READWRITE(F[0]) CONSTREFERENCE(A, A2, F, MMParH, Fi)),
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-R1, N-N2, R1, Fi.mOne, A + R1*lda, lda, A2, lda, Fi.one, F, lda, MMParH));

                    // G <- C2 - E V1
                    TASK(MODE(READ(R1, A[R1], A3[0], MMParH) READWRITE(G[0]) CONSTREFERENCE(Fi, A, A3, G, MMParH)),
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-R1, R1, Fi.mOne, A3, lda, A+R1, lda, Fi.one, G, lda, MMParH));

                    CHECK_DEPENDENCIES;

                    P2 = FFLAS::fflas_new<size_t>(M2-R1);
                    Q2 = FFLAS::fflas_new<size_t>(N-N2);
                    //typename Field::Element * A4R2 = 0;
                    // F = P2 [ L2 ] [ U2 V2 ] Q2
                    //        [ M2 ]
                    TASK(MODE(CONSTREFERENCE(Fi, P2, Q2, F,/* A4R2,*/ R2, PLUQParH1) WRITE(R2/*, A4R2[0]*/) READWRITE(F[0], P2, Q2) ),
                         R2 = PLUQ( Fi, Diag, M2-R1, N-N2, F, lda, P2, Q2, PLUQParH1)
                         );

                    P3 = FFLAS::fflas_new<size_t>(M-M2);
                    Q3 = FFLAS::fflas_new<size_t>(N2-R1);
                    // G = P3 [ L3 ] [ U3 V3 ] Q3
                    //        [ M3 ]
                    TASK(MODE(CONSTREFERENCE(Fi, G, Q3, P3, R3,PLUQParH2) WRITE(R3, P3, Q3) READWRITE(G[0])),
                         R3 = PLUQ( Fi, Diag, M-M2, N2-R1, G, lda, P3, Q3, PLUQParH2));

                    // H <- A4 - ED
                    TASK(MODE(CONSTREFERENCE(Fi, A3, A2, A4, MMParH) READ(M2, N2, R1, A3[0], A2[0]) READWRITE(A4[0])),
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N-N2, R1, Fi.mOne, A3, lda, A2, lda, Fi.one, A4, lda, MMParH));

                    CHECK_DEPENDENCIES;

                    // [ H1 H2 ] <- P3^T H Q2^T
                    // [ H3 H4 ]
                    PermParH.set_numthreads(std::max(nt,1));
                    TASK(MODE(READ(P3, Q2) CONSTREFERENCE(Fi, A4, Q2, P3) READWRITE(A4[0])),
                         applyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, 0, N-N2, A4, lda, Q2, PermParH);
                         applyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, 0, M-M2, A4, lda, P3, PermParH););

                    CHECK_DEPENDENCIES;
                    // [ E1 ] <- P3^T E
                    // [ E2 ]
                    PermParH.set_numthreads(std::max(nt/4,1)); // TODO adjust the value to the size of each block
                    TASK(MODE(READ(P3) CONSTREFERENCE(Fi, P3, A3, PermParH) READWRITE(A3[0])),
                         applyP( Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M-M2, A3, lda, P3, PermParH));

                    // [ M11 ] <- P2^T M1
                    // [ M12 ]
                    TASK(MODE(READ(P2) CONSTREFERENCE(P2, A, Fi, PermParH) READWRITE(A[R1*lda])),
                         applyP(Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2, PermParH));
                    //applyP(Fi, FflasLeft, FflasNoTrans, R1, 0, M2-R1, A+R1*lda, lda, P2);

                    // [ D1 D2 ] <- D Q2^T
                    TASK(MODE(READ(Q2) CONSTREFERENCE(Fi, Q2, A2, PermParH) READWRITE(A2[0])),
                         applyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N-N2, A2, lda, Q2, PermParH));

                    // [ V1 V2 ] <- V1 Q3^T
                    TASK(MODE(READ(Q3) CONSTREFERENCE(Fi, Q3, A, PermParH) READWRITE(A[R1])),
                         applyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, R1, 0, N2-R1, A+R1, lda, Q3, PermParH));

                    //    CHECK_DEPENDENCIES;

                    // I <- H1 U2^-1
                    // K <- H3 U2^-1
                    TASK(MODE(READ(R2, F[0]) CONSTREFERENCE(Fi, A4, F, TRSMParH, R2) READWRITE(A4[0])),
                         ftrsm( Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, R2, Fi.one, F, lda, A4, lda, TRSMParH));

                    CHECK_DEPENDENCIES;

                    typename Field::Element_ptr temp = 0;

                    TASK(MODE(READ(A4[0], R3) READWRITE(temp[0], R2) CONSTREFERENCE(Fi, A4, temp, R2, R3)),
                         temp = FFLAS::fflas_new (Fi, R3, R2);
                         FFLAS::fassign (Fi, R3, R2, A4, lda, temp, R2);
                        );
                    CHECK_DEPENDENCIES;

                    // J <- L3^-1 I (in a temp)
                    TASK(MODE(READ(R2, R3, G[0]) CONSTREFERENCE(Fi, G, temp, R2, R3, TRSMParH) READWRITE(temp[0])),
                         ftrsm( Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, R2, Fi.one, G, lda, temp, R2, TRSMParH););

                    // N <- L3^-1 H2
                    TASK(MODE(READ(R3, R2, G[0]) CONSTREFERENCE(Fi, G, A4, R3, R2, TRSMParH) READWRITE(A4[R2])),
                         ftrsm(Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, R3, N-N2-R2, Fi.one, G, lda, A4+R2, lda, TRSMParH));

                    CHECK_DEPENDENCIES;

                    // O <- N - J V2
                    TASK(MODE(READ(R2, F[R2]) CONSTREFERENCE(Fi, R2, A4, R3, temp, MMParH) READWRITE(A4[R2], temp[0])),
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, R3, N-N2-R2, R2, Fi.mOne, temp, R2, F+R2, lda, Fi.one, A4+R2, lda, MMParH);
                         FFLAS::fflas_delete (temp);
                         temp=0;
                        );

                    typename Field::Element_ptr R = 0;
                    // R <- H4 - K V2
                    TASK(MODE(READ(R2, R3, M2, N2, A4[R3*lda], F[R2]) CONSTREFERENCE(Fi, R, F, R2, R3, MMParH) READWRITE(R[0])),
                         R = A4 + R2 + R3*lda;
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R2, Fi.mOne, A4+R3*lda, lda, F+R2, lda, Fi.one, R, lda, MMParH)
                        );

                    CHECK_DEPENDENCIES;

                    // R <- R - M3 O
                    TASK(MODE(READ(R3, R2, A4[R2], G[R3*lda]) CONSTREFERENCE(Fi, A4, R, R3, R2, G, MMParH) READWRITE(R[0])),
                         fgemm( Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2-R3, N-N2-R2, R3, Fi.mOne, G+R3*lda, lda, A4+R2, lda, Fi.one, R, lda, MMFullParH));

                    CHECK_DEPENDENCIES;

                    /*
                       size_t * P4 = FFLAS::fflas_new<size_t>(M-M2-R3);
                       size_t * Q4 = FFLAS::fflas_new<size_t>(N-N2-R2);
                       */

                    // H4 = P4 [ L4 ] [ U4 V4 ] Q4
                    //         [ M4 ]
                    //TASK(READ(Fi), NOWRITE(R4), READWRITE(R, P4, Q4), PPLUQ, R4, Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4);
                    TASK(MODE(CONSTREFERENCE(Fi, R4, R, P4, Q4, R2, R3, M2, N2,PSHelper) READWRITE(R[0]) WRITE(R4, P4[0], Q4[0])),
                         P4 = FFLAS::fflas_new<size_t>(M-M2-R3);
                         Q4 = FFLAS::fflas_new<size_t>(N-N2-R2);
                         R4 = PLUQ (Fi, Diag, M-M2-R3, N-N2-R2, R, lda, P4, Q4, PSHelper);
                        );
                    CHECK_DEPENDENCIES;

                    PermParH.set_numthreads(std::max(nt/2,1));
                    // [ E21 M31 0 K1 ] <- P4^T [ E2 M3 0 K ]
                    // [ E22 M32 0 K2 ]
                    TASK(MODE(READ(P4[0], R2, R3, M2) CONSTREFERENCE(Fi, P4, A3, R2, R3) READWRITE(A3[R3*lda])),
                         applyP(Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4, PermParH));
                    //applyP( Fi, FflasLeft, FflasNoTrans, N2+R2, 0, M-M2-R3, A3+R3*lda, lda, P4);

                    // [ D21 D22 ]     [ D2 ]
                    // [ V21 V22 ]  <- [ V2 ] Q4^T
                    // [  0   0  ]     [  0 ]
                    // [ O1   O2 ]     [  O ]
                    TASK(MODE(READ(Q4[0], R2, N2, M2, R3) CONSTREFERENCE(Fi, Q4, A2, R2, R3) READWRITE(A2[R2])),
                         applyP( Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4, PermParH));
                    //applyP( Fi, FflasRight, FflasTrans, M2+R3, 0, N-N2-R2, A2+R2, lda, Q4);

                    // P <- Diag (P1 [ I_R1    ] , P3 [ I_R3    ])
                    //               [      P2 ]      [      P4 ]
                    WAIT;
                    //    TASK(MODE(CONSTREFERENCE(P1, P2, P3, P4, R1, R3, MathP, M2) READ(P1, P2, R1, R3, P3, P4, M2) READWRITE(MathP)),
                    MathP = FFLAS::fflas_new<size_t>(M);
                    composePermutationsLLM (MathP, P1, P2, R1, M2);
                    composePermutationsLLM (MathP+M2, P3, P4, R3, M-M2);
                    for (size_t i=M2; i<M; ++i)
                        MathP[i] += M2;
                    /*	 if (R1+R2 < M2)
                         PermApplyS( MathP, 1,1, M2, R1, R2, R3, R4);*/
                    //	 );

                    //CHECK_DEPENDENCIES;

                    //    WAIT;
                    if (R1+R2 < M2){
                        // P <- P S
                        TASK(MODE(CONSTREFERENCE(R1, R2, R3, R4, MathP, M2) READ(R1, R2, R3, R4, M2) READWRITE(MathP[0])),
                             PermApplyS( MathP, 1,1, M2, R1, R2, R3, R4);
                            );

                        // A <-  S^T A
                        PermParH.set_numthreads(std::max(nt,1));
                        TASK(MODE(READ(R1, R2, R3, R4) CONSTREFERENCE(Fi, A, R1, R2, R3, R4) READWRITE(A[0])),
                             MatrixApplyS( Fi, A, lda, N, M2, R1, R2, R3, R4, PermParH) );
                        //MatrixApplyS(Fi, A, lda, N, M2, R1, R2, R3, R4);
                    }

                    // Q<- Diag ( [ I_R1    ] Q1,  [ I_R2    ] Q2 )
                    //            [      Q3 ]      [      P4 ]
                    MathQ = FFLAS::fflas_new<size_t>(N);
                    TASK(MODE(CONSTREFERENCE(Q1, Q2, Q3, Q4, R1, R2) READ(Q1[0], Q2[0], Q3[0], Q4[0], R1, R2) READWRITE(MathQ[0])),
                         composePermutationsLLM (MathQ, Q1, Q3, R1, N2);
                         composePermutationsLLM (MathQ+N2, Q2, Q4, R2, N-N2);
                         for (size_t i=N2; i<N; ++i)
                         MathQ[i] += N2;
                        );
                    CHECK_DEPENDENCIES;

                    if (R1 < N2){
                        // Q <- T Q
                        TASK(MODE(CONSTREFERENCE(R1, R2, R3, R4) READ(R1, R2, R3, R4) READWRITE(MathQ[0])),
                             PermApplyT (MathQ, 1,1,N2, R1, R2, R3, R4););

                        // A <-   A T^T
                        TASK(MODE(READ(R1, R2, R3, R4) CONSTREFERENCE(Fi, A, R1, R2, R3, R4) READWRITE(A[0])),
                             MatrixApplyT(Fi, A, lda, M, N2, R1, R2, R3, R4, PermParH));
                        //			  MatrixApplyT(Fi, A, lda, M, N2, R1, R2, R3, R4);
                    }
                    CHECK_DEPENDENCIES;
                    TASK(MODE(CONSTREFERENCE(MathP, MathQ) READ(MathP[0], MathQ[0]) READWRITE(P[0], Q[0])),
                         MathPerm2LAPACKPerm (Q, MathQ, N);
                         MathPerm2LAPACKPerm (P, MathP, M);
                        );
                    );
                    FFLAS::fflas_delete( MathQ);
                    FFLAS::fflas_delete( MathP);
                    FFLAS::fflas_delete( P1);
                    FFLAS::fflas_delete( P2);
                    FFLAS::fflas_delete( P3);
                    FFLAS::fflas_delete( P4);
                    FFLAS::fflas_delete( Q1);
                    FFLAS::fflas_delete( Q2);
                    FFLAS::fflas_delete( Q3);
                    FFLAS::fflas_delete( Q4);

                    //);



                    return R1+R2+R3+R4;
                    //#endif
    }

    template<class Field>
    inline size_t
    pPLUQ (const Field& Fi, const FFLAS::FFLAS_DIAG Diag,
          size_t M, size_t N,
          typename Field::Element_ptr A, size_t lda, size_t*P, size_t *Q)
    {
        size_t r;
        FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,
                                    FFLAS::StrategyParameter::Threads> PSHelper;
        PAR_BLOCK{
            PSHelper.set_numthreads(NUM_THREADS);
            r = FFPACK::PLUQ (Fi,Diag,M,N,A,lda,P,Q,PSHelper);
        }
        return r;
    }

}// namespace FFPACK

//#endif // OPENMP
#endif // __FFLASFFPACK_ffpack_ppluq_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
