/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2013 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *<ziad.sultan@imag.fr>
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



namespace FFLAS
{


    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element*
    pfgemm(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element * C, const size_t ldc,
           MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Block, StrategyParameter::Threads> > & H){
        {
            H.parseq.set_numthreads( std::min(H.parseq.numthreads(), std::max((size_t)1,(size_t)(m*n/(__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD)))) );

            MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Sequential> SeqH (H);
            size_t sa = (ta==FFLAS::FflasNoTrans)?lda:1;
            size_t sb = (tb==FFLAS::FflasNoTrans)?1:ldb;
            SYNCH_GROUP({FORBLOCK2D(iter,m,n,H.parseq,
                                    TASK( MODE(
                                               READ(A[iter.ibegin()*sa],B[iter.jbegin()*sb])
                                               CONSTREFERENCE(F, SeqH)
                                               READWRITE(C[iter.ibegin()*ldc+iter.jbegin()])),
                                          fgemm( F, ta, tb, iter.iend()-iter.ibegin(), iter.jend()-iter.jbegin(), k, alpha, A+iter.ibegin()*sa, lda, B+iter.jbegin()*sb, ldb, beta, C+iter.ibegin()*ldc+iter.jbegin(), ldc, SeqH););
                                   );
                        });
        }
        return C;


    }

    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element*
    pfgemm(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr AA, const size_t lda,
           const typename Field::ConstElement_ptr BB, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element * C, const size_t ldc,
           MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::ThreeDAdaptive> > & H){

        typename Field::Element a = alpha;
        typename Field::Element b = beta;
        typename Field::ConstElement_ptr B = BB;
        typename Field::ConstElement_ptr A = AA;
        if (!m || !n) {return C;}
        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }

        if (H.parseq.numthreads()<=1 || std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
            MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
            return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);
        }

        typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::ThreeDAdaptive> > MMH_t;
        MMH_t H1(H);
        MMH_t H2(H);
        if(__FFLASFFPACK_DIMKPENALTY*m > k && m >= n) {
            SYNCH_GROUP(size_t M2= m>>1;
                        H1.parseq.set_numthreads(H1.parseq.numthreads() >> 1);
                        H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());

                        typename Field::ConstElement_ptr A1= A;
                        typename Field::ConstElement_ptr A2= A+M2*((ta==FFLAS::FflasTrans)?1:lda);
                        typename Field::Element_ptr C1= C;
                        typename Field::Element_ptr C2= C+M2*ldc;

                        // 2 multiply (1 split on dimension m)

                        TASK(MODE(CONSTREFERENCE(F, H1) READ(A1,B) READWRITE(C1)),
                             {pfgemm( F, ta, tb, M2, n, k, alpha, A1, lda, B, ldb, beta, C1, ldc, H1);}
                            );

                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,B) READWRITE(C2)),
                             {pfgemm(F, ta, tb, m-M2, n, k, alpha, A2, lda, B, ldb, beta, C2, ldc, H2);}
                            );
                       );

        } else if (__FFLASFFPACK_DIMKPENALTY*n > k) {
            SYNCH_GROUP(
                        size_t N2 = n>>1;
                        H1.parseq.set_numthreads( H1.parseq.numthreads() >> 1);
                        H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());
                        typename Field::ConstElement_ptr B1= B;
                        typename Field::ConstElement_ptr B2= B+N2*((tb==FFLAS::FflasTrans)?ldb:1);

                        typename Field::Element_ptr C1= C;
                        typename Field::Element_ptr C2= C+N2;

                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A,B1) READWRITE(C1)), pfgemm(F, ta, tb, m, N2, k, a, A, lda, B1, ldb, b, C1, ldc, H1));
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A,B2) READWRITE(C2)), pfgemm(F, ta, tb, m, n-N2, k, a, A, lda, B2, ldb, b,C2, ldc, H2));
                       );

        } else {

            size_t K2 = k>>1;

            typename Field::ConstElement_ptr B1= B;
            typename Field::ConstElement_ptr B2= B+K2*((tb==FFLAS::FflasTrans)?1:ldb);
            typename Field::ConstElement_ptr A1= A;
            typename Field::ConstElement_ptr A2= A+K2*((ta==FFLAS::FflasTrans)?lda:1);
            typename Field::Element_ptr C2 = fflas_new (F, m, n,Alignment::CACHE_PAGESIZE);

            H1.parseq.set_numthreads(H1.parseq.numthreads() >> 1);
            H2.parseq.set_numthreads(H.parseq.numthreads()-H1.parseq.numthreads());
            SYNCH_GROUP(
                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A1,B1) READWRITE(C)), pfgemm(F, ta, tb, m, n, K2, a, A1, lda, B1, ldb, b, C, ldc, H1));

                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,B2) READWRITE(C2)), pfgemm(F, ta, tb, m, n, k-K2, a, A2, lda, B2, ldb, F.zero, C2, n, H2));
                        CHECK_DEPENDENCIES;

                        TASK(MODE(CONSTREFERENCE(F) READ(C2) READWRITE(C)),faddin(F, n, m, C2, n, C, ldc));

                       );
            fflas_delete(C2);
        }

        return C;
    }

    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element*
    pfgemm (const Field& F,
            const FFLAS_TRANSPOSE ta,
            const FFLAS_TRANSPOSE tb,
            const size_t m,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            const typename Field::ConstElement_ptr AA, const size_t lda,
            const typename Field::ConstElement_ptr BB, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element * C, const size_t ldc,
            MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::TwoDAdaptive> > & H){

        typename Field::Element a = alpha;
        typename Field::Element b = beta;
        typename Field::ConstElement_ptr B = BB;
        typename Field::ConstElement_ptr A = AA;
        if (!m || !n) {return C;}
        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }
        if (H.parseq.numthreads()<=1 || m*n<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
            MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
            return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);

        }
        typedef MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::TwoDAdaptive> > MMH_t;
        MMH_t H1(H);
        MMH_t H2(H);
        H1.parseq.set_numthreads(H1.parseq.numthreads() >> 1);
        H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());
        if(m >= n) {
            size_t M2= m>>1;
            typename Field::ConstElement_ptr A1= A;
            typename Field::ConstElement_ptr A2= A+M2*((ta==FFLAS::FflasTrans)?1:lda);
            typename Field::Element_ptr C1= C;
            typename Field::Element_ptr C2= C+M2*ldc;
            SYNCH_GROUP(
                        TASK(MODE(CONSTREFERENCE(F,H1, A1, B) READ(M2, A1[0],B[0]) READWRITE(C1[0])), pfgemm(F, ta, tb, M2, n, k, alpha, A1, lda, B, ldb, beta, C1, ldc, H1));
                        TASK(MODE(CONSTREFERENCE(F,H2, A2, B) READ(M2, A2[0],B[0]) READWRITE(C2[0])), pfgemm(F, ta, tb, m-M2, n, k, alpha, A2, lda, B, ldb, beta, C2, ldc, H2));

                       );

        } else {
            size_t N2 = n>>1;
            typename Field::ConstElement_ptr B1= B;
            typename Field::ConstElement_ptr B2= B+N2*((tb==FFLAS::FflasTrans)?ldb:1);
            typename Field::Element_ptr C1= C;
            typename Field::Element_ptr C2= C+N2;
            SYNCH_GROUP(
                        TASK(MODE(CONSTREFERENCE(F,H1, A, B1) READ(N2, A[0], B1[0]) READWRITE(C1[0])), pfgemm(F, ta, tb, m, N2, k, a, A, lda, B1, ldb, b, C1, ldc, H1));
                        TASK(MODE(CONSTREFERENCE(F,H2, A, B2) READ(N2, A[0], B2[0]) READWRITE(C2[0])), pfgemm(F, ta, tb, m, n-N2, k, a, A, lda, B2, ldb, b,C2, ldc, H2));
                       );
        }
        return C;
    }

    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element*
    pfgemm( const Field& F,
            const FFLAS_TRANSPOSE ta,
            const FFLAS_TRANSPOSE tb,
            const size_t m,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            const typename Field::ConstElement_ptr AA, const size_t lda,
            const typename Field::ConstElement_ptr BB, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element * C, const size_t ldc,
            MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::TwoD> > & H){

        typename Field::Element a = alpha;
        typename Field::Element b = beta;
        typename Field::ConstElement_ptr B = BB;
        typename Field::ConstElement_ptr A = AA;
        if (!m || !n) {return C;}
        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }

        if(H.parseq.numthreads()<=1|| m*n<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
            MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Sequential> SeqH(H);
            return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, SeqH);
        } else
        {
            size_t M2= m>>1;
            size_t N2= n>>1;

            typename Field::ConstElement_ptr A1= A;
            typename Field::ConstElement_ptr A2= A+M2*((ta==FFLAS::FflasTrans)?1:lda);
            typename Field::ConstElement_ptr B1= B;
            typename Field::ConstElement_ptr B2= B+N2*((tb==FFLAS::FflasTrans)?ldb:1);

            typename Field::Element_ptr C11= C;
            typename Field::Element_ptr C21= C+M2*ldc;
            typename Field::Element_ptr C12= C+N2;
            typename Field::Element_ptr C22= C+N2+M2*ldc;

            typedef MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::TwoD> > MMH_t;
            MMH_t H1(H);
            MMH_t H2(H);
            MMH_t H3(H);
            MMH_t H4(H);
            size_t nt = H.parseq.numthreads();
            size_t nt_rec = nt/4;
            size_t nt_mod = nt%4;
            H1.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H2.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H3.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H4.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            SYNCH_GROUP(
                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A1,B1) READWRITE(C11)), pfgemm(F, ta, tb, M2, N2, k, alpha, A1, lda, B1, ldb, beta, C11, ldc, H1));

                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A1,B2) READWRITE(C12)), pfgemm(F, ta, tb, M2, n-N2, k, alpha, A1, lda, B2, ldb, beta, C12, ldc, H2));

                        TASK(MODE(CONSTREFERENCE(F,H3) READ(A2,B1) READWRITE(C21)), pfgemm(F, ta, tb, m-M2, N2, k, a, A2, lda, B1, ldb, b, C21, ldc, H3));

                        TASK(MODE(CONSTREFERENCE(F,H4) READ(A2,B2) READWRITE(C22)), pfgemm(F, ta, tb, m-M2, n-N2, k, a, A2, lda, B2, ldb, b,C22, ldc, H4));
                       );
        }
        return C;
    }



    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element_ptr
    pfgemm(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m,
           const size_t n,
           const size_t k,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr B, const size_t ldb,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc,
           MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::ThreeD> > & H){


        if (!m || !n) {return C;}
        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }
        if(H.parseq.numthreads() <= 1|| std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){
            FFLAS::MMHelper<Field, AlgoT, FieldTrait,FFLAS::ParSeqHelper::Sequential> WH (H);
            return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, WH);
        }
        else
        {
            typename Field::Element a = alpha;
            typename Field::Element b = 0;

            size_t M2= m>>1;
            size_t N2= n>>1;
            size_t K2= k>>1;
            typename Field::ConstElement_ptr A11= A;
            typename Field::ConstElement_ptr A12= A+K2*((ta==FFLAS::FflasTrans)?lda:1);
            typename Field::ConstElement_ptr A21= A+M2*((ta==FFLAS::FflasTrans)?1:lda);
            typename Field::ConstElement_ptr A22= A12+M2*((ta==FFLAS::FflasTrans)?1:lda);

            typename Field::ConstElement_ptr B11= B;
            typename Field::ConstElement_ptr B12= B+N2*((tb==FFLAS::FflasTrans)?ldb:1);
            typename Field::ConstElement_ptr B21= B+K2*((tb==FFLAS::FflasTrans)?1:ldb);
            typename Field::ConstElement_ptr B22= B12+K2*((tb==FFLAS::FflasTrans)?1:ldb);

            typename Field::Element_ptr C11= C;
            typename Field::Element_ptr C_11 = fflas_new (F, M2, N2,Alignment::CACHE_PAGESIZE);

            typename Field::Element_ptr C12= C+N2;
            typename Field::Element_ptr C_12 = fflas_new (F, M2, n-N2,Alignment::CACHE_PAGESIZE);

            typename Field::Element_ptr C21= C+M2*ldc;
            typename Field::Element_ptr C_21 = fflas_new (F, m-M2, N2,Alignment::CACHE_PAGESIZE);

            typename Field::Element_ptr C22= C+N2+M2*ldc;
            typename Field::Element_ptr C_22 = fflas_new (F, m-M2, n-N2,Alignment::CACHE_PAGESIZE);

            // 1/ 8 multiply in parallel
            //omp_set_task_affinity(omp_get_locality_domain_num_for( C11));

            typedef MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::ThreeD> > MMH_t;
            MMH_t H1(H);
            MMH_t H2(H);
            MMH_t H3(H);
            MMH_t H4(H);
            MMH_t H5(H);
            MMH_t H6(H);
            MMH_t H7(H);
            MMH_t H8(H);
            size_t nt = H.parseq.numthreads();
            size_t nt_rec = nt/8;
            size_t nt_mod = nt % 8 ;
            H1.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H2.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H3.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H4.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H5.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H6.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H7.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H8.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));

            SYNCH_GROUP(
                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A11,B11) READWRITE(C11)), pfgemm(F, ta, tb, M2, N2, K2, alpha, A11, lda, B11, ldb, beta, C11, ldc, H1));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C_11));
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A12,B21) WRITE(C_11)), pfgemm(F, ta, tb, M2, N2, k-K2, a, A12, lda, B21, ldb, b,C_11, N2, H2));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C12));
                        TASK(MODE(CONSTREFERENCE(F,H3) READ(A12,B22) READWRITE(C12)), pfgemm(F, ta, tb, M2, n-N2, k-K2, alpha, A12, lda, B22, ldb, beta, C12, ldc, H3));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C_12));
                        TASK(MODE(CONSTREFERENCE(F,H4) READ(A11,B12) WRITE(C_12)), pfgemm(F, ta, tb, M2, n-N2, K2, a, A11, lda, B12, ldb, b, C_12, n-N2, H4));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C21));
                        TASK(MODE(CONSTREFERENCE(F,H5) READ(A22,B21) READWRITE(C21)), pfgemm(F, ta, tb, m-M2, N2, k-K2, alpha, A22, lda, B21, ldb, beta, C21, ldc, H5));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C_21));
                        TASK(MODE(CONSTREFERENCE(F,H6) READ(A21,B11) WRITE(C_21)), pfgemm(F, ta, tb, m-M2, N2, K2, a, A21, lda, B11, ldb, b,C_21, N2, H6));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C22));
                        TASK(MODE(CONSTREFERENCE(F,H7) READ(A21,B12) READWRITE(C22)), pfgemm(F, ta, tb, m-M2, n-N2, K2, alpha, A21, lda, B12, ldb, beta, C22, ldc, H7));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C_22));
                        TASK(MODE(CONSTREFERENCE(F,H8) READ(A22,B22) WRITE(C_22)), pfgemm(F, ta, tb, m-M2, n-N2, k-K2, a, A22, lda, B22, ldb, b,C_22, n-N2, H8));

                        CHECK_DEPENDENCIES;
                        // 2/ final add
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C11));
                        TASK(MODE(CONSTREFERENCE(F) READ(C_11) READWRITE(C11)), faddin(F, M2, N2, C_11, N2, C11, ldc));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C12));
                        TASK(MODE(CONSTREFERENCE(F) READ(C_12) READWRITE(C12)),faddin(F, M2, n-N2, C_12, n-N2, C12, ldc));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C21));
                        TASK(MODE(CONSTREFERENCE(F) READ(C_21) READWRITE(C21)), faddin(F, m-M2, N2, C_21, N2, C21, ldc));
                        //omp_set_task_affinity(omp_get_locality_domain_num_for( C22));
                        TASK(MODE(CONSTREFERENCE(F) READ(C_22) READWRITE(C22)), faddin(F, m-M2, n-N2, C_22, n-N2, C22, ldc));

                        );
                        FFLAS::fflas_delete (C_11);
                        FFLAS::fflas_delete (C_12);
                        FFLAS::fflas_delete (C_21);
                        FFLAS::fflas_delete (C_22);
        }
        return C;
    }

    template<class Field, class AlgoT, class FieldTrait>
    typename Field::Element*
    pfgemm( const Field& F,
            const FFLAS_TRANSPOSE ta,
            const FFLAS_TRANSPOSE tb,
            const size_t m,
            const size_t n,
            const size_t k,
            const typename Field::Element alpha,
            const typename Field::ConstElement_ptr A, const size_t lda,
            const typename Field::ConstElement_ptr B, const size_t ldb,
            const typename Field::Element beta,
            typename Field::Element_ptr C, const size_t ldc,
            MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::ThreeDInPlace> > & H){


        if (!m || !n) {return C;}
        if (!k || F.isZero (alpha)){
            fscalin(F, m, n, beta, C, ldc);
            return C;
        }

        if(H.parseq.numthreads() <= 1|| std::min(m*n,std::min(m*k,k*n))<=__FFLASFFPACK_SEQPARTHRESHOLD*__FFLASFFPACK_SEQPARTHRESHOLD){	// threshold
            FFLAS::MMHelper<Field, AlgoT, FieldTrait,FFLAS::ParSeqHelper::Sequential> WH (H);
            return fgemm(F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, WH);
        }else{
            size_t M2= m>>1;
            size_t N2= n>>1;
            size_t K2= k>>1;
            typename Field::ConstElement_ptr A11= A;
            typename Field::ConstElement_ptr A12= A+K2*((ta==FFLAS::FflasTrans)?lda:1);
            typename Field::ConstElement_ptr A21= A+M2*((ta==FFLAS::FflasTrans)?1:lda);
            typename Field::ConstElement_ptr A22= A12+M2*((ta==FFLAS::FflasTrans)?1:lda);

            typename Field::ConstElement_ptr B11= B;
            typename Field::ConstElement_ptr B12= B+N2*((tb==FFLAS::FflasTrans)?ldb:1);
            typename Field::ConstElement_ptr B21= B+K2*((tb==FFLAS::FflasTrans)?1:ldb);
            typename Field::ConstElement_ptr B22= B12+K2*((tb==FFLAS::FflasTrans)?1:ldb);


            typename Field::Element_ptr C11= C;
            typename Field::Element_ptr C12= C+N2;
            typename Field::Element_ptr C21= C+M2*ldc;
            typename Field::Element_ptr C22= C+N2+M2*ldc;
            typedef MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive,StrategyParameter::ThreeDInPlace> > MMH_t;
            MMH_t H1(H);
            MMH_t H2(H);
            MMH_t H3(H);
            MMH_t H4(H);
            size_t nt = H.parseq.numthreads();
            size_t nt_rec = nt/4;
            size_t nt_mod = nt%4;
            H1.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H2.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H3.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            H4.parseq.set_numthreads(std::max(size_t(1),nt_rec + ((nt_mod-- > 0)?1:0)));
            SYNCH_GROUP(
                        // 1/ 4 multiply
                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A11,B11) READWRITE(C11)), pfgemm(F, ta, tb, M2, N2, K2, alpha, A11, lda, B11, ldb, beta, C11, ldc, H1));
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A12,B22) READWRITE(C12)), pfgemm(F, ta, tb, M2, n-N2, k-K2, alpha, A12, lda, B22, ldb, beta, C12, ldc, H2));
                        TASK(MODE(CONSTREFERENCE(F,H3) READ(A22,B21) READWRITE(C21)), pfgemm(F, ta, tb, m-M2, N2, k-K2, alpha, A22, lda, B21, ldb, beta, C21, ldc, H3));
                        TASK(MODE(CONSTREFERENCE(F,H4) READ(A21,B12) READWRITE(C22)), pfgemm(F, ta, tb, m-M2, n-N2, K2, alpha, A21, lda, B12, ldb, beta, C22, ldc, H4));

                        CHECK_DEPENDENCIES;
                        // 2/ 4 add+multiply
                        TASK(MODE(CONSTREFERENCE(F,H1) READ(A12,B21) READWRITE(C11)), pfgemm(F, ta, tb, M2, N2, k-K2, alpha, A12, lda, B21, ldb, F.one, C11, ldc, H1));
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A11,B12) READWRITE(C12)), pfgemm(F, ta, tb, M2, n-N2, K2, alpha, A11, lda, B12, ldb, F.one, C12, ldc, H2));
                        TASK(MODE(CONSTREFERENCE(F,H3) READ(A21,B11) READWRITE(C21)), pfgemm(F, ta, tb, m-M2, N2, K2, alpha, A21, lda, B11, ldb, F.one, C21, ldc, H3));
                        TASK(MODE(CONSTREFERENCE(F,H4) READ(A22,B22) READWRITE(C22)), pfgemm(F, ta, tb, m-M2, n-N2, k-K2, alpha, A22, lda, B22, ldb, F.one, C22, ldc, H4));
                       );
        }
        return C;
    }



} // FFLAS
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
