/* fflas/fflas_pftrsm.inl
 * Copyright (C) 2013 Ziad Sultan
 *
 * Written by Ziad Sultan  < Ziad.Sultan@imag.fr >
 * Time-stamp: <18 Dec 15 16:09:24 Jean-Guillaume.Dumas@imag.fr>
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


#ifndef __FFLASFFPACK_fflas_pftrsm_INL
#define __FFLASFFPACK_fflas_pftrsm_INL

#define PTRSM_HYBRID_THRESHOLD 256

#include "fflas-ffpack/paladin/parallel.h"

namespace FFLAS {

    template<class Field, class Cut, class Param>
    inline typename Field::Element_ptr
    ftrsm( const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const FFLAS::FFLAS_UPLO UpLo,
           const FFLAS::FFLAS_TRANSPOSE TA,
           const FFLAS::FFLAS_DIAG Diag,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           TRSMHelper <StructureHelper::Iterative, ParSeqHelper::Parallel<Cut,Param> > & H)
    // const FFLAS::CuttingStrategy method,
    // const size_t numThreads)
    {
        typedef TRSMHelper<StructureHelper::Recursive,ParSeqHelper::Sequential> seqRecHelper;
        SYNCH_GROUP(
                    seqRecHelper SeqH(H);
                    if(Side == FflasRight){
                    FORBLOCK1D(iter, m, H.parseq,
                               TASK(MODE(READ(A[0]) CONSTREFERENCE(F, A, B, SeqH,H) READWRITE(B[iter.begin()*ldb])), ftrsm( F, Side, UpLo, TA, Diag, iter.end()-iter.begin(), n, alpha, A, lda, B + iter.begin()*ldb, ldb, SeqH));
                              );
                    } else {
                    FORBLOCK1D(iter, n, H.parseq,
                               //							  seqRecHelper SeqH(H);
                               TASK(MODE(READ(A[0]) CONSTREFERENCE(F, A, B, SeqH,H) READWRITE(B[iter.begin()])), ftrsm(F, Side, UpLo, TA, Diag, m, iter.end()-iter.begin(), alpha, A , lda, B + iter.begin(), ldb, SeqH));
                              );
                    }
                   );
        return B;
    }
    template<class Field, class Cut, class Param>
    inline typename Field::Element_ptr
    ftrsm( const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const FFLAS::FFLAS_UPLO UpLo,
           const FFLAS::FFLAS_TRANSPOSE TA,
           const FFLAS::FFLAS_DIAG Diag,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
#ifdef __FFLAS__TRSM_READONLY
           typename Field::ConstElement_ptr
#else
           typename Field::Element_ptr
#endif
           A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           TRSMHelper <StructureHelper::Hybrid, ParSeqHelper::Parallel<Cut, Param> > & H)
    // const FFLAS::CuttingStrategy method,
    // const size_t numThreads)
    {

        if(Side == FflasRight){

            size_t nt = H.parseq.numthreads();
            size_t nt_it,nt_rec;
            if (m/PTRSM_HYBRID_THRESHOLD < nt){
                nt_it = (int)ceil(double(m)/PTRSM_HYBRID_THRESHOLD);
                nt_rec = (int)ceil(double(nt)/double(nt_it));
            } else { nt_it = nt; nt_rec = 1;}
            //			ForStrategy1D<size_t> iter(m, ParSeqHelper::Parallel((size_t)nt_it,H.parseq.method));
            //			for (iter.begin(); ! iter.end(); ++iter) {
            //			SYNCH_GROUP(H.parseq.numthreads(),
            SYNCH_GROUP(
                        ParSeqHelper::Parallel<Cut,Param> psh(nt_rec);
                        TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Parallel<Cut,Param> > SeqH (psh);
                        H.parseq.set_numthreads(nt_it);
                        FORBLOCK1D(iter, m, H.parseq,
                                   //				      std::cerr<<"trsm_rec nt = "<<nt_rec<<std::endl;
                                   TASK(MODE(READ(A) CONSTREFERENCE(F, A, B, SeqH) READWRITE(B[iter.begin()*ldb])),
                                        ftrsm( F, Side, UpLo, TA, Diag, iter.end()-iter.begin(), n, alpha, A, lda, B + iter.begin()*ldb, ldb, SeqH));
                                  );
                       );

        } else {

            size_t nt = H.parseq.numthreads();
            size_t nt_it=nt;
            size_t nt_rec=1;
            if (nt_it*PTRSM_HYBRID_THRESHOLD >= n){
                nt_it>>=1;
                nt_rec<<=1;
                while(nt_it*PTRSM_HYBRID_THRESHOLD >= n){
                    nt_it>>=1;
                    nt_rec<<=1;
                }
                nt_it<<=1;
                nt_rec>>=1;
            }

            // if ((int)n/PTRSM_HYBRID_THRESHOLD < nt){
            // 	nt_it = std::min(nt,(int)ceil(double(n)/PTRSM_HYBRID_THRESHOLD));
            // 	nt_rec = ceil(double(nt)/nt_it);
            // } else { nt_it = nt; nt_rec = 1;}

            //	ForStrategy1D<size_t> iter(n, ParSeqHelper::Parallel((size_t)nt_it,H.parseq.method));
            //				for (iter.begin(); ! iter.end(); ++iter) {


            // 				std::cerr<<"trsm_rec nt_it = "<<nt_it<<std::endl;
            // 				std::cerr<<"trsm_rec nt_rec = "<<nt_rec<<std::endl;

            SYNCH_GROUP(
                        ParSeqHelper::Parallel<Cut,Param> psh(nt_rec);
                        TRSMHelper<StructureHelper::Recursive, ParSeqHelper::Parallel<Cut,Param> > SeqH (psh);
                        H.parseq.set_numthreads(nt_it);
                        FORBLOCK1D(iter, n, H.parseq,
                                   TASK(MODE(READ(A) CONSTREFERENCE(F, A, B, SeqH) READWRITE(B[iter.begin()])), ftrsm( F, Side, UpLo, TA, Diag, m, iter.end()-iter.begin(), alpha, A , lda, B + iter.begin(), ldb, SeqH));
                                  );
                       );
        }
        return B;
        }

        } // FFLAS


#endif // __FFLASFFPACK_fflas_pftrsm_INL
        /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
        // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
