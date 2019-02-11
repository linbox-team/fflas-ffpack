/* fflas-ffpack/ffpack/ffpack_charpoly_kgfast.inl
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

#ifndef __FFLASFFPACK_ffpack_charpoly_kgfastgeneralized_INL
#define __FFLASFFPACK_ffpack_charpoly_kgfastgeneralized_INL


//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using
// Keller-Gehrig's fast algorithm.
//---------------------------------------------------------------------


#include <iostream> // std::cout
#include "fflas-ffpack/utils/fflas_io.h"
#ifdef __FFLASFFPACK_DEBUG
namespace FFPACK {
    template <class Field>
    typename Field::Element_ptr buildMatrix (const Field& F,
                                             typename Field::ConstElement_ptr E,
                                             typename Field::ConstElement_ptr C,
                                             const size_t lda,
                                             const size_t*B,
                                             const size_t*T,
                                             const size_t me,
                                             const size_t mc,
                                             const size_t lambda,
                                             const size_t mu);
    template <class Field>
    void printA(const Field& F,
                std::ostream& os,
                typename Field::ConstElement_ptr E,
                typename Field::ConstElement_ptr C,
                const size_t lda,
                const size_t*B,
                const size_t*T,
                const size_t me,const size_t mc, const size_t lambda, const size_t mu)
    {

        typename Field::Element_ptr A = buildMatrix(F,E,C,lda,B,T,me,mc,lambda,mu);
        size_t N = mc+me+lambda+mu;
        FFLAS::WriteMatrix(os,F,N,N,A,N);
        FFLAS::fflas_delete (A);
    }
} // FFPACK
#endif

namespace FFPACK {
    template <class Field>
    typename Field::Element_ptr buildMatrix (const Field& F,
                                             typename Field::ConstElement_ptr E,
                                             typename Field::ConstElement_ptr C,
                                             const size_t lda,
                                             const size_t*B,
                                             const size_t*T,
                                             const size_t me,
                                             const size_t mc,
                                             const size_t lambda,
                                             const size_t mu)
    {

        size_t N = mc+me+lambda+mu;
        typename Field::Element_ptr A = FFLAS::fflas_new (F, N, N);
        for (size_t j=0; j<lambda+me;++j)
            if (B[j] < N){
                for (size_t i=0;i<N;++i)
                    F.assign( *(A+i*N+j), F.zero);
                F.assign( *(A+B[j]*lda+j), F.one);
            } else {
                FFLAS::fassign (F, N, E+B[j]-N, lda, A+j, N);
            }
        for (size_t j=lambda+me; j<lambda+me+mu; ++j)
            for (size_t i=0;i<N;++i)
                F.assign( *(A+i*N+j), F.zero);
        for (size_t i=0; i<mu; ++i)
            F.assign( *(A+(lambda+me+mc+i)*lda+lambda+me+T[i]), F.one);
        for (size_t j=0; j<mc; ++j)
            FFLAS::fassign(F,N,C+j,lda,A+N-mc+j,N);
        //! @bug is this :
        // FFLAS::fassign(F,N,mc,C,lda,A+N-mc,N);
        return A;
    }

    namespace Protected {

        template <class Field, class Polynomial>
        std::list<Polynomial>&
        KGFast_generalized (const Field& F, std::list<Polynomial>& charp,
                            const size_t N,
                            typename Field::Element_ptr A, const size_t lda)
        {

            //std::cerr<<"Dans KGFast"<<std::endl;
            size_t mc=N>>1; // Matrix A is transformed into a mc_Frobenius form
            size_t me=N-mc;
            // B[i] = j, the row of the 1 if the col Ai is sparse;
            // B[i] = n+k, if the col Ai is the kth col of E
            size_t * B = FFLAS::fflas_new<size_t>(N);
            bool * allowedRows = FFLAS::fflas_new<bool>(N);
            for (size_t i=0;i<(N+1)/2;++i)
                allowedRows[i]=true;
            // T[i] = j si T_i,j = 1
            size_t * T = FFLAS::fflas_new<size_t>(N);
            for (size_t i=0;i<N;++i)
                T[i]=i;
            size_t lambda=0;

            typename Field::Element_ptr C, E = A;
#ifdef __FFLASFFPACK_DEBUG
            std::cerr<<"Debut KGFG"<<std::endl
            <<" ----------------------------"<<std::endl;
#endif
            int exit_value = 0 ;
            while (mc > 0) {
#ifdef __FFLASFFPACK_DEBUG
                std::cerr<<"Boucle1: mc,me,lambda="<<mc<<" "<<me<<" "<<lambda<<std::endl;
                //FFLAS::WriteMatrix (std::cerr, F, N, N, A, lda);
#endif
                size_t mu=0;
                C = A + (N-mc);
                for (size_t i = 0; i<me;++i)
                    B[lambda+i] = N+i;
#ifdef __FFLASFFPACK_DEBUG
                for (size_t i=0;i<lambda+me;++i)
                    std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
                //std::cerr<<std::endl<<"mc="<<mc<<":";
#endif
                while (mu < N-mc && !exit_value) {
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"Boucle2: mu,me,lambda="<<mu<<" "<<me<<" "<<lambda<<std::endl;
                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // B1 <- C1^-1.B1
                    std::cerr<<"Forming LUP";
#endif
                    size_t ncols = ((mu==0)||(mc<=mu))?mc:mc-mu;
                    typename Field::Element_ptr LUP = FFLAS::fflas_new (F, lambda+me, ncols);
                    for (size_t i=0;i < lambda + me; ++i)
                        if (allowedRows[i])
                            FFLAS::fassign (F, ncols, C+i*lda, 1, LUP+i*ncols, 1);
                        else
                            for (size_t j = 0; j < ncols; ++j)
                                F.assign (*(LUP+i*ncols+j), F.zero);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
                    //FFLAS::WriteMatrix (std::cerr<<"LUP="<<std::endl,F,lambda+me,ncols,LUP,ncols);
                    std::cerr<<"LQUP(C1)";
#endif
                    size_t * P = FFLAS::fflas_new<size_t>(ncols);
                    size_t * Q = FFLAS::fflas_new<size_t>(lambda+me);
                    for (size_t i=0; i<ncols;++i)
                        P[i]=0;
                    for (size_t i=0; i<lambda+me;++i)
                        Q[i]=0;

                    size_t r = LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, lambda + me, ncols, LUP, ncols, P, Q);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
#endif

                    if (r==0){
                        if ((lambda == 0) && (ncols == mc)){
                            std::cerr<<"BLOCAGE lambda=0!!!"<<std::endl;
                            //Rec call on the leading block
                            KGFast_generalized (F, charp, me, A, lda);

                            //Rec call on the trailing block
                            typename Field::Element_ptr At = buildMatrix(F,E,C,lda,B,T,me,mc,lambda,mu);
                            KGFast_generalized (F, charp, N-me, At+me*(lda+1), lda);
                            FFLAS::fflas_delete (At);
                            exit_value = -1;
                            break;

                        } else if (me != 0) {
                            std::cerr<<"BLOCAGE me!=0!!!"<<std::endl;
                            exit_value = -1;
                            break ;

                        }
                        else {
                            for (int i = int(mu+1); i--; )
                                T[(size_t)i+lambda] = T[i]+lambda;
                            for (size_t i=0; i< lambda; ++i)
                                T[B[i]-mc-1] = i;
                            mu += lambda;
                            lambda = 0;
                            break;
                        }
                        //std::cerr<<"BLOCAGE !!!"<<std::endl;
                        //exit(-1);
                    }

#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"Forming genreric rank profil C1";
                    // form the generic rank profil block C1 Q^TPAP^TQ
                    for (size_t i=0;i<r;++i)
                        std::cerr<<"P["<<i<<"] = "<<P[i]<<std::endl;
#endif
                    applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
                            N, 0, (int)r, C, lda, P);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
#endif
                    //printA(F,cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
                    // (E, C) <- P(E, C)
                    applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                            me, 0, (int)r, E+(N-mc)*lda, lda, P);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
                    //printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif
                    applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                            mc, 0, (int)r, C+(N-mc)*lda, lda, P);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
#endif
                    //printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
                    // T <- P T

                    // !!!!!!!Attention -> ajouter le traitement du cas 0<mu<mc
                    for (size_t k = 0; k<r; ++k)
                        if (P[k] > (size_t) k){
                            if ((mu>=mc-k)){
#ifdef __FFLASFFPACK_DEBUG
                                std::cerr<<"// on permute LN-mc+k et L_N-mc+P[k]"<<std::endl;
#endif
                                size_t tmp = T[mu-mc+k];
                                T[mu-mc+k] = T[mu-mc+P[k]];
                                T[mu-mc+P[k]] = tmp;
                            }
                            else if (mu){
                                std::cerr<<"CAS MU < MC - k"<<std::endl;
                                exit_value = -1;
                                break;
                            }
                            // Updating B to be improved (tabulated B^-1)
                            for (size_t i=0; i<lambda+me; ++i){
                                if (B[i] == N-mc+k)
                                    B[i] = N-mc+P[k];
                                else if (B[i] == N-mc+P[k])
                                    B[i] = N-mc+k;
                            }

                        }
                    if (exit_value)
                        break;
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
                    //printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif

                    // (E, C) <- Q^T(E, C)
                    applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                            me, 0,(int) r, E, lda, Q);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
                    //printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif
                    applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                            mc, 0, (int)r, C, lda, Q);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
#endif
                    // F <- Q^T F
                    size_t * tempP = FFLAS::fflas_new<size_t>(lambda+me+mc);
                    for (size_t i=0; i< lambda+me+mc; ++i)
                        tempP[i] = i;

                    for (int i = int(r) ; i--; )
                        if (Q[i] > (size_t) i){
#ifdef __FFLASFFPACK_DEBUG
                            std::cerr<<"Permutation de tempP["<<i
                            <<"] et tempP["<<Q[i]<<"]"<<std::endl;
#endif
                            // on permute LN-mc+k et L_N-mc+P[k]
                            size_t tmp = tempP[i];
                            tempP[i] = tempP[Q[i]];
                            tempP[Q[i]] = tmp;
                        }

#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
#endif
                    for (size_t i=0; i < lambda+me; ++i)
                        if (B[i] < N)
                            B[i] = tempP[B[i]];
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<".";
#endif
                    FFLAS::fflas_delete( tempP);

#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<std::endl<<"Avant B<-BQ"<<std::endl;
                    for (size_t i=0; i<lambda+me;++i)
                        std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
#endif
                    // B <- B Q
                    for (int k = int(r-1); k>=0; --k)
                        if (Q[k] > (size_t) k){
                            // on permute Ck et C_Q[k]
                            size_t tmp = B[k];
                            B[k] = B[Q[k]];
                            B[Q[k]] = tmp;
                        }
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"Apres"<<std::endl;
                    for (size_t i=0; i<lambda+me;++i)
                        std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;

                    std::cerr<<".";
#endif

                    // grouping the bloc L in LUP
                    for (size_t i=0; i<r; ++i)
                        if (Q[i]>i)
                            FFLAS::fassign(F, i, LUP+Q[i]*mc,1, LUP+i*mc, 1);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

#if 0
                    std::cerr<<"LUP="<<std::endl;
                    //FFLAS::WriteMatrix (std::cerr, F, mc, mc, LUP, mc);
                    std::cerr<<" "<<r;
#endif
                    // E'1 <- C11^-1 E1
                    std::cerr<<"// E'1 <- C11^-1 E1";
#endif

                    ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                          r, me, F.one, LUP, mc , E, lda);
                    ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                          r, me, F.one, LUP, mc , E, lda);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // C'12 <- C11^-1 C12
                    std::cerr<<"// C'12 <- C11^-1 C12";
#endif
                    ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit,
                          r, mc-r, F.one, LUP, mc , C+r, lda);
                    ftrsm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
                          r, mc-r, F.one, LUP, mc , C+r, lda);
                    FFLAS::fflas_delete (LUP);
                    FFLAS::fflas_delete( P);
                    FFLAS::fflas_delete( Q);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif

                    // E'2 <- E2 - C21.E'1
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"// E'2 <- E2 - C21.E'1";
#endif
                    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-r, me, r,
                          F.mOne, C+r*lda, lda, E, lda,
                          F.one, E+r*lda, lda);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
                    // C'22 <- C22 - C21.C'12
                    std::cerr<<"// C'22 <- C22 - C21.C'12";
#endif
                    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-r, mc-r, r,
                          F.mOne, C+r*lda, lda, C+r, lda,
                          F.one, C+r*(lda+1), lda);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // Shifting E: E1;E2 -> E2;E1
                    std::cerr<<"// Shifting E: E1;E2 -> E2;E1";
#endif
                    typename Field::Element_ptr tmp = FFLAS::fflas_new (F, r, me);
                    for (size_t i=0; i<r; ++i)
                        FFLAS::fassign (F, me, E+i*lda, 1, tmp+i*me, 1);
                    for (size_t i=r; i< N; ++i)
                        FFLAS::fassign (F, me, E+i*lda, 1, E+(i-r)*lda, 1);
                    for (size_t i=0; i<r; ++i)
                        FFLAS::fassign (F, me, tmp+i*me, 1, E+(i+N-r)*lda, 1);
                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    // Shifting C_{*,2}: C_{1,2};C_{2,2} -> C_{2,2};C_{1,2}
                    std::cerr<<"// Shifting C_{*,2}: C_{1,2};C_{2,2} -> C_{2,2};C_{1,2}";
#endif
                    tmp = FFLAS::fflas_new (F, r, mc-r);
                    for (size_t i=0; i<r; ++i)
                        FFLAS::fassign (F, mc-r, C+r+i*lda, 1, tmp+i*(mc-r), 1);
                    for (size_t i=r; i< N; ++i)
                        FFLAS::fassign (F, mc-r, C+r+i*lda, 1, C+r+(i-r)*lda, 1);
                    for (size_t i=0; i<r; ++i)
                        FFLAS::fassign (F, mc-r, tmp+i*(mc-r), 1, C+r+(i+N-r)*lda, 1);
                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // C'2 <- T C2
                    std::cerr<<"// C'2 <- T C2";
#endif
                    // To be improved!!!
                    tmp = FFLAS::fflas_new (F, mu, r);
                    typename Field::Element_ptr C2 = C+(N-mu-mc)*lda;
                    for (size_t i=0; i<mu; ++i)
                        FFLAS::fassign (F, r, C2+T[i]*lda, 1, tmp+i*r, 1);
                    for (size_t i=0; i<mu; ++i)
                        FFLAS::fassign (F, r, tmp+i*r, 1, C2+i*lda, 1);
                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    // [C'2;C'3] += [E2;E3].C
                    std::cerr<<"// [C'2;C'3] += [E2;E3].C";
#endif
                    tmp = FFLAS::fflas_new (F, me, r);
                    for (size_t i=0; i<lambda+me; ++i)
                        if (B[i] >= N){
                            FFLAS::fassign (F, r, C+i*lda, 1, tmp+(B[i]-N)*r, 1);
                        }
                    fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mu + r, r, me,
                           F.one, E+(N-mu-r)*lda, lda, tmp, r,
                           F.one, C+(N-mu-mc)*lda, lda);

                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    // shifting [C'2;C'3]
                    std::cerr<<"// shifting [C'2;C'3]";
#endif
                    tmp = FFLAS::fflas_new (F, mc-r, r);
                    typename Field::Element_ptr C4 = C + (N-mc+r)*lda;
                    for (size_t i=0; i < (mc-r); ++i){
                        FFLAS::fassign (F, r, C4 + i*lda, 1, tmp+i*r, 1);
                    }
                    for (int i = int(N-1); i >= (int) (N -mu-r); --i)
                        FFLAS::fassign (F, r, C+((size_t)i-mc+r)*lda, 1, C+i*(int)lda, 1);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);


                    // tmp2 <- C'1 (the rows corresponding to E)
                    std::cerr<<"// tmp2 <- C'1 (the rows corresponding to E)";
#endif
                    typename Field::Element_ptr tmp2 = FFLAS::fflas_new (F, me, r);
                    for (size_t i = 0; i < lambda+me; ++i)
                        if (B[i] >= N){
#ifdef __FFLASFFPACK_DEBUG
                            std::cerr<<"saving in row "<<B[i]-N<<std::endl;
#endif
                            FFLAS::fassign (F, r, C+i*lda, 1, tmp2+(B[i]-N)*r, 1);
                        }
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    // C'_F[i] <- C_i
                    std::cerr<<"// C'_F[i] <- C_i";
                    std::cerr<<"lambda,r,me = "<<lambda<<" "<<r<<" "<<me<<std::endl;
#endif
                    typename Field::Element_ptr tmp3 = FFLAS::fflas_new (F, lambda+me,r);

                    for (size_t i = 0; i < lambda+me; ++i)
                        if (B[i] < N){
#ifdef __FFLASFFPACK_DEBUG
                            std::cerr<<"copie de la ligne "<<i<<std::endl;
#endif
                            FFLAS::fassign (F, r, C + i*lda, 1, tmp3 + i*r, 1);
                        }
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"1"<<std::endl;
#endif
                    for (size_t i = 0; i < N-mu-r; ++i)
                        for (size_t j = 0; j < r; ++j)
                            F.assign (*(C+i*lda+j), F.zero);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"2"<<std::endl;
#endif
                    for (size_t i = 0; i < lambda+me; ++i){
#ifdef __FFLASFFPACK_DEBUG
                        std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
#endif
                        if (B[i] < N)
                            FFLAS::fassign (F, r, tmp3+i*r, 1, C+(B[i]-r)*lda, 1);
                    }
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"3"<<std::endl;
#endif
                    FFLAS::fflas_delete (tmp3);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // C'1 += E1 tmp2
                    std::cerr<<"// C'1 += E1 tmp2";
#endif
                    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N-mu-r, r, me,
                          F.one, E, lda, tmp2, r, F.one, C, lda);
                    FFLAS::fflas_delete (tmp2);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // C'_1 += C_2 C4
                    std::cerr<<"// C'_1 += C_2 C4";
#endif
                    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N, r, mc-r,
                          F.one, C+r, lda, tmp, r, F.one, C, lda);
                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;

                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

                    // switching C_1 <-> C_2
                    std::cerr<<"// switching C_1 <-> C_2";
#endif
                    tmp = FFLAS::fflas_new (F, N, r);
                    for (size_t j = 0; j<r; ++j)
                        FFLAS::fassign (F, N, C+j, lda, tmp+j, r);
                    for (size_t j = r; j<mc; ++j)
                        FFLAS::fassign (F, N, C+j, lda, C+j-r, lda);
                    for (size_t j = 0; j<r; ++j)
                        FFLAS::fassign (F, N, tmp+j, r, C+mc-r+j, lda);
                    FFLAS::fflas_delete (tmp);
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;


                    printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);


                    // update the datastructure:
                    std::cerr<<"// update the datastructure:";
#endif
                    mu += r;
                    tmp2 = FFLAS::fflas_new (F, N, me);
                    size_t nlambda= 0, nme=0;
                    for (size_t i=0;i<lambda+me;++i)
                        allowedRows[i]=true;
                    for (size_t j=r; j < lambda + me; ++j){
                        if (B[j] >= N){
#ifdef __FFLASFFPACK_DEBUG
                            std::cerr<<"B["<<j-r<<"] = "<<N+nme<<std::endl;
#endif
                            FFLAS::fassign (F, N, E+(B[j]-N), lda, tmp2+nme, me);
                            B[j-r] = N + nme;
                            nme++;
                        } else {
#ifdef __FFLASFFPACK_DEBUG
                            std::cerr<<"B["<<j-r<<"] = "<<B[j]<<std::endl;
#endif
                            B[j-r] = B[j]-r;
                            allowedRows[B[j]-r] = false;
                            nlambda++;
                        }
                    }
                    for (size_t j=0; j<nme; ++j)
                        FFLAS::fassign (F, N, tmp2+j, me, E+j, lda);
                    lambda = nlambda;
                    me = nme;
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"..done"<<std::endl;
#endif
                    FFLAS::fflas_delete (tmp2);
                }
                // update the datastructure: F <- T
                for (size_t i=0; i<mu; ++i){
#ifdef __FFLASFFPACK_DEBUG
                    std::cerr<<"B[T["<<i<<"]] = "<<"B["<<T[i]<<"] = "<<mc+i<<std::endl;
#endif

                    B[T[i]] = mc+i;
                    T[i]=i;
                }
                E=C;
                me = mc;
                mc>>=1;
                me -= mc;
                lambda = mu;
                for (size_t i=0;i<me+mc;++i)
                    allowedRows[i]=true;
                for (size_t i=me+mc;i<lambda+me+mc;++i)
                    allowedRows[i]=false;

            }

            FFLAS::fflas_delete( B );
            FFLAS::fflas_delete( T );
            FFLAS::fflas_delete( allowedRows );

            if (exit_value)
                exit(exit_value);
            Polynomial *minP = new Polynomial();
            minP->resize(N+1);
            minP->operator[](N) = F.one;
            typename Polynomial::iterator it = minP->begin();
            for (size_t j=0; j<N; ++j, it++){
                F.neg(*it, *(A+N-1+j*lda));
            }
            charp.clear();
            charp.push_front(*minP);
            return charp;
        }
    } // Protected
} // FFPACK

#endif // __FFLASFFPACK_ffpack_charpoly_kgfastgeneralized_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
