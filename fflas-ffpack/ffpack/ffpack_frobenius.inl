/* fflas-ffpack/ffpack/ffpack_frobenius.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
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

#include <givaro/givranditer.h>

//---------------------------------------------------------------------
// CharpolyArithProg: Las Vegas algorithm to compute the Charpoly
// over a large field (Z/pZ, s.t.  p > 2n^2)
//---------------------------------------------------------------------
//
//
namespace FFPACK { namespace Protected {
    template <class Field>
    inline void CompressRows (Field& F, const size_t M,
                       typename Field::Element_ptr A, const size_t lda,
                       typename Field::Element_ptr tmp, const size_t ldtmp,
                       const size_t * d, const size_t nb_blocs);

    template <class Field>
    inline void CompressRowsQK (Field& F, const size_t M,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d,const size_t deg, const size_t nb_blocs);

    template <class Field>
    inline void DeCompressRows (Field& F, const size_t M, const size_t N,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d, const size_t nb_blocs);

    template <class Field>
    inline void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr tmp, const size_t ldtmp,
                           const size_t * d, const size_t deg, const size_t nb_blocs);

    template <class Field>
    inline void CompressRowsQA (Field& F, const size_t M,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d, const size_t nb_blocs);

    template <class Field>
    inline void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr tmp, const size_t ldtmp,
                           const size_t * d, const size_t nb_blocs);


    template <class PolRing>
    inline void
    RandomKrylovPrecond (const PolRing& PR, std::list<typename PolRing::Element>& completedFactors, const size_t N,
                         typename PolRing::Domain_t::Element_ptr A, const size_t lda,
                         size_t & Nb, typename PolRing::Domain_t::Element_ptr& B, size_t& ldb,
                         typename PolRing::Domain_t::RandIter& g, const size_t degree)
    {
        typedef typename PolRing::Domain_t Field;
        typedef typename PolRing::Element Polynomial;
        const Field& F = PR.getdomain();

        FFLASFFPACK_check(degree);
        size_t noc = (N-1)/degree + 1;
        // Building the workplace matrix
        typename Field::Element_ptr K  = FFLAS::fflas_new (F, degree*noc, N);
        typename Field::Element_ptr K2 = FFLAS::fflas_new (F, degree*noc, N);
        size_t ldk = N;
//        size_t bk_stride = noc*ldk;

        size_t *dA = FFLAS::fflas_new<size_t>(N); //PA
        size_t *dK = FFLAS::fflas_new<size_t>(noc*degree);
        for (size_t i=0; i<noc; ++i)
            dK[i]=0;
#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        Givaro::Timer timkrylov, timelim, timsimilar, timrest;
        timkrylov.start();
#endif
            // Picking a random noc x N block vector U^T
        Givaro::GeneralRingNonZeroRandIter<Field> nzg (g);
        RandomMatrix (F, noc, N, K, ldk*degree, g);
        for (size_t i = 0; i < noc; ++i)
            nzg.random (*(K + i*(degree*ldk+1)));

        // Computing the bloc Krylov matrix [u1 Au1 .. A^(c-1) u1 u2 Au2 ...]^T
        for (size_t i = 1; i<degree; ++i){
            fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans,  noc, N, N,F.one,
                   K+(i-1)*ldk, degree*ldk, A, lda, F.zero, K+i*ldk, degree*ldk);
        }
        // K2 <- K (re-ordering)
        //! @todo swap to save space ??
        //! @todo don't assing K2 c*noc x N but only mas (c,noc) x N and store each one after the other 
        // size_t w_idx = 0;
        // for (size_t i=0; i<noc; ++i)
        //     for (size_t j=0; j<degree; ++j, w_idx++)
        //         FFLAS::fassign(F, N, (K+(i+j*noc)*ldk), 1, (K2+(w_idx)*ldk), 1);

        // Copying K2 <- K
//        for (size_t i=0; i<noc*degree; ++i)
        FFLAS::fassign (F, noc*degree, N, K, ldk, K2, ldk);

#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        timkrylov.stop();
         std::cerr<<"Random Krylov Preconditionning:"<<std::endl
                 <<"  Krylov basis computation : "<<timkrylov.usertime()<<std::endl;
        timelim.start();
#endif
        size_t * Pk = FFLAS::fflas_new<size_t>(N);
        size_t nrows = noc*degree;
        size_t * Qk = FFLAS::fflas_new<size_t>(nrows);
        for (size_t i=0; i<nrows; ++i)
            Qk[i] = 0;
        for (size_t i=0; i<N; ++i)
            Pk[i] = 0;

            // @todo: replace by PLUQ
        size_t R = LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, noc*degree, N, K, ldk, Pk, Qk);
        size_t row_idx = 0;
        size_t ii=0;
        size_t dold = degree;
        size_t nb_full_blocks = 0;
        size_t Mk = 0;
        // Determining the degree sequence dK
        for (size_t k = 0; k<noc; ++k){
            size_t d = 0;
            while ( (d<degree) && (row_idx<R) && (Qk[row_idx] == ii)) {ii++; row_idx++; d++;}
            if (d > dold){
                // std::cerr << "FAIL in preconditionning phase:"
                //           << " degree sequence is not monotonically not increasing"
                // 	     << std::endl;
                FFLAS::fflas_delete (K, Pk, Qk, dA, dK);
                throw CharpolyFailed();
            }
            dK[k] = dold = d;
            Mk++;
            if (d == degree)
                nb_full_blocks++;
            if (row_idx < N)
                ii = Qk[row_idx];
        }
#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        timelim.stop();
        std::cerr <<"  LU (Krylov)              : "<<timelim.usertime()<<std::endl;
        timsimilar.start();
#endif
        // Selection of the last iterate of each block

        typename Field::Element_ptr K3 = FFLAS::fflas_new (F, Mk, N);
        typename Field::Element_ptr K4 = FFLAS::fflas_new (F, Mk, N);
        size_t bk_idx = 0;
        for (size_t i = 0; i < Mk; ++i){
            FFLAS::fassign (F, N, (K2 + (bk_idx + dK[i]-1)*ldk), 1, (K3+i*ldk), 1);
            bk_idx += degree;
        }
        FFLAS::fflas_delete (K2);

        // K <- K A^T
        fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans, Mk, N, N,F.one,  K3, ldk, A, lda, F.zero, K4, ldk);

        // K <- K P^T
        applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans, Mk, 0,(int) R, K4, ldk, Pk);

        // K <- K U^-1
        ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, Mk, R,F.one, K, ldk, K4, ldk);

        // L <-  Q^T L
        applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N, 0,(int) R, K, ldk, Qk);

        // K <- K L^-1
        ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Mk, R,F.one, K, ldk, K4, ldk);

        //undoing permutation on L
        applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans, N, 0,(int) R, K, ldk, Qk);

        // Recovery of the completed invariant factors
        size_t Ma = Mk;
        size_t Ncurr = R;
        size_t offset = Ncurr-1;
        // How many non-full blocks are completed :
        //  - if no full_blocks: all of them
        //  - if there is at least one full block, then the first non-full will be further eliminated later on and should not be pushed as completed
        size_t last_completed = nb_full_blocks + nb_full_blocks?1:0;
        for (size_t ip1=Mk; ip1 > last_completed;  --ip1){
                size_t i = ip1-1;
                if (dK[i] >= 1){
                for (size_t j = offset+1; j<R; ++j)
                    if (!F.isZero(*(K4 + i*ldk + j))){
                        //std::cerr<<"FAIL C != 0 in preconditionning"<<std::endl;
                        FFLAS::fflas_delete (K,K3,K4,Pk,Qk,dA,dK);
                        throw CharpolyFailed();
                    }
                Polynomial P (dK [i]+1);
                F.assign(P[dK[i]],F.one);
                for (size_t j=0; j < dK [i]; ++j)
                    F.neg (P [dK [i]-j-1], *(K4 + i*ldk + (offset-j)));
                completedFactors.push_front(P);
                offset -= dK [i];
                Ncurr -= dK [i];
                Ma--;
            }
        }
        Mk = Ma;

#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        timsimilar.stop();
        std::cerr <<"  Similarity)              : "<<timsimilar.usertime()<<std::endl;
        timrest.start();
#endif
        if (R<N){
                // The Krylov basis did not span the whole space
                // Recurse on the complementary subspace
            if (! FFLAS::fiszero (F, nb_full_blocks+1, N-R, K4+R, ldk)){
                
            // for (size_t i=0; i<nb_full_blocks + 1; ++i)
            //     for (size_t j=R; j<N; ++j){
            //         if (!F.isZero( *(K4+i*ldk+j) )){
                        FFLAS::fflas_delete (K3, K4, K, Pk, Qk, dA, dK);
                        throw CharpolyFailed();
            }

            size_t Nrest = N-R;
            typename Field::Element_ptr K21 = K + R*ldk;
            typename Field::Element_ptr K22 = K21 + R;

            //  Compute the n-k last rows of A' = P A^T P^T in K2_
            // A = A . P^t
            applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
                    N, 0,(int) R, A, lda, Pk);

            // Copy K2_ = (A'_2)^t
            for (size_t i=0; i<Nrest; i++)
                FFLAS::fassign (F, N, A+R+i, lda, K21+i*ldk, 1);
            
            // A = A . P : Undo the permutation on A
            applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, N, 0,(int) R, A, lda, Pk);

            // K2_ = K2_ . P^t (=  ( P A^t P^t )2_ )
            applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, Nrest, 0,(int) R, K21, ldk, Pk);

            // K21 = K21 . S1^-1
            ftrsm (F, FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit, Nrest, R, F.one, K, ldk, K21, ldk);

            typename Field::Element_ptr Arec = FFLAS::fflas_new (F, Nrest, Nrest);
            size_t ldarec = Nrest;

            // Creation of the matrix A2 for recursive call
            FFLAS::fassign (F, Nrest, Nrest, K22, ldk, Arec, ldarec);

            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, R,F.mOne, K21, ldk, K+R, ldk,F.one, Arec, ldarec);

            std::list<Polynomial> polyList;
            polyList.clear();

            // Recursive call on the complementary subspace
            CharPoly (PR, polyList, Nrest, Arec, ldarec, g, FfpackArithProgKrylovPrecond);
            FFLAS::fflas_delete (Arec);
            completedFactors.merge(polyList);
        }
#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
            timrest.stop();
            std::cerr<<"  left-over                : "<<timrest.usertime()<<std::endl;
#endif

        FFLAS::fflas_delete (K, K3, Pk, Qk);
        for (size_t i=0; i<Mk; ++i)
            dA[i] = dK[i];
        bk_idx = 0;

        ldb = Ma;
        Nb = Ncurr;
        B = FFLAS::fflas_new (F, Ncurr, ldb);
        for (size_t j=0; j<Ma; ++j)
            FFLAS::fassign(F, Ncurr, K4+j*ldk, 1, B+j, ldb);
        FFLAS::fflas_delete (dA, dK, K4);

    }

    template <class PolRing>
    inline std::list<typename PolRing::Element>&
    ArithProg (const PolRing& PR, std::list<typename PolRing::Element>& frobeniusForm,
               const size_t N, typename PolRing::Domain_t::Element_ptr A, const size_t lda,
               const size_t degree)
    {
        if (!N) return frobeniusForm;
        typedef typename PolRing::Domain_t Field;
        typedef typename PolRing::Element Polynomial;
        const Field& F = PR.getdomain();

#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        Givaro::Timer tim;
        tim.start();
#endif
	size_t nb_full_blocks, Mk, Ma;
	Mk = Ma = nb_full_blocks = (N-1)/degree +1;
	typename Field::Element_ptr K, K3, Ac;
	size_t ldk=Ma;
	size_t ldac = Ma;
	K = FFLAS::fflas_new(F, N, ldk);
	Ac = FFLAS::fflas_new(F, N, ldac);

	FFLAS::fassign(F, N, Ma, A, lda, Ac, ldac);
	FFLAS::fassign(F, N, Ma, Ac, ldac, K, ldk);

        size_t Ncurr=N;

        size_t * dA = FFLAS::fflas_new<size_t>(Ma);
	size_t * dK = FFLAS::fflas_new<size_t>(Mk);
	for (size_t i=0; i<Ma; i++){
		dK[i] = dA[i] = degree;
	}
	size_t rdeg = N % degree;
	if (rdeg)
		dK[Mk-1] = dA[Ma-1] = rdeg;

        typename Field::Element_ptr Arp = FFLAS::fflas_new (F, Ncurr, Ma);
        size_t ldarp = Ncurr;
        size_t deg = degree+1;
        size_t * rp=NULL;
            // Main loop of the arithmetic progession
        while ((nb_full_blocks >= 1) && (Mk > 1)) {
            size_t block_idx, it_idx, rp_val;
            FFLAS::fflas_delete (K);
            K = FFLAS::fflas_new (F, Ncurr, Ma);
            K3 = FFLAS::fflas_new (F, Ncurr, Ma);
            ldk = Ma;

            // Computation of the rank profile
            for (size_t i=0; i < Ncurr; ++i)
                for (size_t j=0; j < Ma; ++j)
                    *(Arp + j*ldarp + Ncurr-i-1) = *(Ac + i*ldac + j);
            rp = FFLAS::fflas_new<size_t>(2*Ncurr);
            for (size_t i=0; i<2*Ncurr; ++i)
                rp[i] = 0;
            size_t RR;
            try{
                RR = SpecRankProfile (F, Ma, Ncurr, Arp, ldarp, deg-1, rp);
            } catch (CharpolyFailed){
                FFLAS::fflas_delete (Arp, Ac, K, K3, rp, dA, dK);
                throw CharpolyFailed();
            }
            if (RR < Ncurr){
                //std::cerr<<"FAIL RR<Ncurr"<<std::endl;
                FFLAS::fflas_delete (Arp, Ac, K, K3, rp, dA, dK);
                throw CharpolyFailed();
            }

            // Computation of the degree vector dK
            it_idx = 0;
            rp_val = 0;
            size_t gg = 0;
            size_t dtot=0;
            block_idx = 0;
            nb_full_blocks = 0;
            while (dtot<Ncurr){
                do {gg++; rp_val++; it_idx++;}
                while ( /*(gg<Ncurr ) &&*/ (rp[gg] == rp_val) && (it_idx < deg ));
                if ((block_idx)&&(it_idx > dK[block_idx-1])){
                    FFLAS::fflas_delete (Arp, Ac, K, K3, rp, dA, dK);
                    throw CharpolyFailed();
                    //std::cerr<<"FAIL d non decreasing"<<std::endl;
                    //exit(-1);
                }
                dK[block_idx++] = it_idx;
                dtot += it_idx;
                if (it_idx == deg)
                    nb_full_blocks ++;
                it_idx=0;
                rp_val = rp[gg];
            }

            Mk = block_idx;

            // Selection of dense colums of K
            for (size_t i=0; i < nb_full_blocks; ++i){
                FFLAS::fassign (F, Ncurr, Ac+i, ldac, K+i, ldk);
            }

            // K <- QK K
            size_t pos = nb_full_blocks*(deg-1);
            for (size_t i = nb_full_blocks; i < Mk; ++i){
                for (size_t j=0; j<Ncurr; ++j)
                    F.assign (*(K + i + j*ldk), F.zero);
                F.assign (*(K + i + (pos + dK[i]-1)*ldk),F.one);
                pos += dA[i];
            }

            // Copying K3 <- K
            for (size_t i=0; i<Mk; ++i)
                FFLAS::fassign (F, Ncurr, K+i, ldk, K3+i, ldk);
            Protected::CompressRowsQK (F, Mk, K3 + nb_full_blocks*(deg-1)*ldk, ldk,
                                       Arp, ldarp, dK+nb_full_blocks, deg, Mk-nb_full_blocks);

            // K <- PA K
            Protected::CompressRows (F, nb_full_blocks, K, ldk, Arp, ldarp, dA, Ma);

            // A <- newQA^T K (compress)
            Protected::CompressRowsQA (F, Ma, Ac, ldac, Arp, ldarp, dA, Ma);

            // K <- A K
            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ncurr-Ma, nb_full_blocks, Ma,F.one,
                   Ac, ldac, K+(Ncurr-Ma)*ldk, ldk,F.one, K, ldk);
            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ma, nb_full_blocks, Ma,F.one,
                   Ac+(Ncurr-Ma)*ldac, ldac, K+(Ncurr-Ma)*ldk, ldk, F.zero, Arp, ldarp);
            for (size_t i=0; i< Ma; ++i)
                FFLAS::fassign(F, nb_full_blocks, Arp+i*ldarp, 1, K+(Ncurr-Ma+i)*ldk, 1);

            // Copying the last rows of A times K
            size_t offset = (deg-2)*nb_full_blocks;
            for (size_t i = nb_full_blocks; i < Mk; ++i) {
                for (size_t j=0; j<Ncurr; ++j)
                    F.assign(*(K+i+j*ldk), F.zero);
                if (dK[i] == dA[i]) // copy the column of A
                    FFLAS::fassign (F, Ncurr, Ac+i, ldac, K+i, ldk);
                else{
                    F.assign (*(K + i + (offset+dK[i]-1)*ldk),F.one);
                }
                offset += dA[i]-1;
            }

            // K <- QA K
            Protected::DeCompressRowsQA (F, Mk, Ncurr, K, ldk, Arp, ldarp, dA, Ma);

            // K <- QK^T K
            Protected::CompressRowsQK (F, Mk, K + nb_full_blocks*(deg-1)*ldk, ldk, Arp, ldarp,
                                       dK+nb_full_blocks, deg, Mk-nb_full_blocks);

            // K <- K^-1 K
            size_t *P=FFLAS::fflas_new<size_t>(Mk);
            size_t *Q=FFLAS::fflas_new<size_t>(Mk);
            if (LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, Mk, Mk , K3 + (Ncurr-Mk)*ldk, ldk, P, Q) < Mk){
                // should never happen (not a LAS VEGAS check)
                //std::cerr<<"FAIL R2 < MK"<<std::endl;
                //			exit(-1);
            }
            ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Mk, Mk,F.one,
                   K3 + (Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
            ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, Mk, Mk,F.one,
                   K3+(Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
            applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
                    Mk, 0,(int) Mk, K+(Ncurr-Mk)*ldk,ldk, P);
            fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ncurr-Mk, Mk, Mk,F.mOne,
                   K3, ldk, K+(Ncurr-Mk)*ldk,ldk,F.one, K, ldk);
            FFLAS::fflas_delete( P);
            FFLAS::fflas_delete( Q);

            // K <- PK^T K
            Protected::DeCompressRows (F, Mk, Ncurr, K, ldk, Arp, ldarp, dK, Mk);

            // K <- K PK (dA <- dK)
            if (nb_full_blocks*deg < Ncurr)
                Ma = nb_full_blocks+1;
            else
                Ma = nb_full_blocks;

            for (size_t i=0; i< Ma; ++i)
                dA[i] = dK[i];

            // Recovery of the completed invariant factors
            offset = Ncurr-1;
            size_t oldNcurr = Ncurr;
            for (size_t i=Mk-1; i>=nb_full_blocks+1;  --i)
                if (dK[i] >= 1){
                    Polynomial  PP (dK [i]+1);
                    F.assign(PP[dK[i]],F.one);
                    for (size_t j=0; j < dK[i]; ++j)
                        F.neg( PP[dK[i]-j-1], *(K + i + (offset-j)*ldk));
                    frobeniusForm.push_front(PP);
                    offset -= dK[i];
                    Ncurr -= dK[i];
                }
            for (size_t i= offset+1; i<oldNcurr; ++i)
                for (size_t j=0; j<nb_full_blocks+1; ++j){
                    if (!F.isZero( *(K+i*ldk+j) )){
                        //std::cerr<<"FAIL C != 0"<<std::endl;
                        FFLAS::fflas_delete (rp, Arp, Ac, K, K3, dA, dK);
                        throw CharpolyFailed();
                    }
                }

            // A <- K
            FFLAS::fflas_delete (Ac);
            FFLAS::fflas_delete (Arp);
            Ac = FFLAS::fflas_new (F, Ncurr, Mk);
            ldac = Mk;
            Arp = FFLAS::fflas_new (F, Ncurr, Mk);
            ldarp=Ncurr;
            for (size_t i=0; i < Ncurr; ++i )
                FFLAS::fassign (F, Mk, K + i*ldk, 1, Ac + i*ldac, 1);

            deg++;
            FFLAS::fflas_delete (K3, rp);
        }

        // Recovery of the first invariant factor
        Polynomial Pl(dK [0]+1);
        F.assign(Pl[dK[0]],F.one);
        for (size_t j=0; j < dK[0]; ++j)
            F.neg( Pl[j], *(K  + j*ldk));
        frobeniusForm.push_front(Pl);
        FFLAS::fflas_delete (Arp, Ac, K, dA, dK);
#ifdef __FFLASFFPACK_ARITHPROG_PROFILING
        tim.stop();
        std::cerr<<"Arith Prog                 : "<<tim.usertime()<<std::endl;
#endif
        return frobeniusForm;
    }

    template <class Field>
    void CompressRowsQK (Field& F, const size_t M,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d, const size_t deg,const size_t nb_blocs)
    {

        int currtmp = 0;
        size_t currw = d[0]-1;
        size_t currr = d[0]-1;
        for (int i = 0; i< int(nb_blocs)-1; ++i){
            // FFLAS::fassign(F,deg-d[i],M,A+currr*lda,lda,tmp+(size_t)currtmp*ldtmp);
            for (int j = int(d[i]-1); j<int(deg)-1; ++j, ++currr, ++currtmp)
                FFLAS::fassign(F, M,  A + currr*lda, 1, tmp + (size_t)currtmp*ldtmp, 1);
            // currr += (deg - d[i]);
            for (int j=0; j < int(d[i+1]) -1; ++j, ++currr, ++currw){
                FFLAS::fassign(F, M, A+(currr)*lda, 1, A + (currw)*lda, 1);
            }
        }
        for (int i=0; i < currtmp; ++i, ++currw){
            FFLAS::fassign (F, M, tmp + (size_t)i*ldtmp, 1, A + (currw)*lda, 1);
        }
    }

    template <class Field>
    void CompressRows (Field& F, const size_t M,
                       typename Field::Element_ptr A, const size_t lda,
                       typename Field::Element_ptr tmp, const size_t ldtmp,
                       const size_t * d, const size_t nb_blocs)
    {

        size_t currd = d[0]-1;
        size_t curri = d[0]-1;
        for (int i = 0; i< int(nb_blocs)-1; ++i){
            FFLAS::fassign(F, M,  A + currd*lda, 1, tmp + i*(int)ldtmp, 1);
            for (int j=0; j < int(d[i+1]) -1; ++j){
                FFLAS::fassign(F, M, A+(currd+(size_t)j+1)*lda, 1, A + (curri++)*lda, 1);
            }
            currd += d[i+1];
        }
        for (int i=0; i < int(nb_blocs)-1; ++i){
            FFLAS::fassign (F, M, tmp + i*(int)ldtmp, 1, A + (curri++)*lda, 1);
        }
    }

    template <class Field>
    void DeCompressRows (Field& F, const size_t M, const size_t N,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d, const size_t nb_blocs)
    {

        for (int i=0; i<int(nb_blocs)-1; ++i)
            FFLAS::fassign(F, M, A + (N-nb_blocs+(size_t)i)*lda, 1, tmp + i*(int)ldtmp, 1);

        size_t w_idx = N - 2;
        size_t r_idx = N - nb_blocs - 1;
        int i = int(nb_blocs)-1 ;
        for (; i--; ){
            for (size_t j = 0; j<d[i+1]-1; ++j)
                FFLAS::fassign (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
            FFLAS::fassign (F, M, tmp + i*(int)ldtmp, 1, A + (w_idx--)*lda, 1);
        }
    }

    template <class Field>
    void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr tmp, const size_t ldtmp,
                           const size_t * d, const size_t deg,const size_t nb_blocs)
    {

        size_t zeroblockdim = 1; // the last block contributes with 1
        size_t currtmp = 0;
        for (int i=0; i<int(nb_blocs)-1; ++i)
            zeroblockdim += deg - d[i];
        for (size_t i=0; i < zeroblockdim - 1; ++i, ++currtmp)
            FFLAS::fassign(F, M,  A + (N - zeroblockdim +i)*lda, 1, tmp + currtmp*ldtmp, 1);
        currtmp--;
        size_t w_idx = N - 2;
        size_t r_idx = N - zeroblockdim - 1;

        int i = int(nb_blocs)-1 ;
        for (; i--;){
            for (size_t j = 0; j < d [i+1] - 1; ++j)
                FFLAS::fassign (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
            for (size_t j = 0; j < deg - d[i]; ++j)
                FFLAS::fassign (F, M, tmp + (currtmp--)*ldtmp, 1, A + (w_idx--)*lda, 1);
        }
    }

    template <class Field>
    void CompressRowsQA (Field& F, const size_t M,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr tmp, const size_t ldtmp,
                         const size_t * d, const size_t nb_blocs)
    {

        size_t currd = 0;
        size_t curri = 0;
        for (size_t i = 0; i< nb_blocs; ++i){
            FFLAS::fassign(F, M,  A + currd*lda, 1, tmp + i*ldtmp, 1);
            for (size_t j=0; j < d[i] -1; ++j)
                FFLAS::fassign(F, M, A+(currd+j+1)*lda, 1, A + (curri++)*lda, 1);
            currd += d[i];
        }
        for (size_t i=0; i < nb_blocs; ++i)
            FFLAS::fassign (F, M, tmp + i*ldtmp, 1, A + (curri++)*lda, 1);
    }

    template <class Field>
    void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr tmp, const size_t ldtmp,
                           const size_t * d, const size_t nb_blocs)
    {

        for (size_t i=0; i<nb_blocs; ++i)
            FFLAS::fassign(F, M, A + (N-nb_blocs+i)*lda, 1, tmp + i*ldtmp, 1);

        size_t w_idx = N - 1;
        size_t r_idx = N - nb_blocs - 1;
        int i = int(nb_blocs) ;
        for (; i--; ){
            for (size_t j = 0; j<d[i]-1; ++j)
                FFLAS::fassign (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
            FFLAS::fassign (F, M, tmp + i*(int)ldtmp, 1, A + (w_idx--)*lda, 1);
        }
    }

} // Protected
} //FFPACK
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
