/* fflas/fflas_pfgemv.inl
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
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > & H){

        if (H.parseq.numthreads()==1 || m <= 1){
            fgemv(F, ta,  m, n,  alpha, A, lda, X, incX, beta, Y, incY);

        }else{
            typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > MMH_t;
            MMH_t H1(H);
            MMH_t H2(H);
            size_t Ydim = (ta == FflasNoTrans) ? m : n;
            size_t Ydim2 = Ydim>>1;
            SYNCH_GROUP(
                        H1.parseq.set_numthreads(H.parseq.numthreads() >> 1);
                        H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());
                        typename Field::ConstElement_ptr A1 = A;
                        typename Field::ConstElement_ptr A2 = A + ((ta==FflasNoTrans) ? Ydim2*lda : Ydim);
                        typename Field::Element_ptr C1 = Y;
                        typename Field::Element_ptr C2 = Y + Ydim2*incY;
                        if (ta==FflasNoTrans){
                        TASK(CONSTREFERENCE(F,H1) MODE( READ(A1,X) READWRITE(C1)),
                             {fgemv( F, ta,  Ydim2, n, alpha, A1, lda, X, incX, beta, C1, incY, H1);}
                            );
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,X) READWRITE(C2)),
                             {fgemv(F, ta, m-Ydim2, n, alpha, A2, lda, X, incX, beta, C2, incY, H2);}
                            );
                        } else {
                        TASK(CONSTREFERENCE(F,H1) MODE( READ(A1,X) READWRITE(C1)),
                             {fgemv( F, ta,  m, Ydim2, alpha, A1, lda, X, incX, beta, C1, incY, H1);}
                            );
                        TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,X) READWRITE(C2)),
                             {fgemv(F, ta, m, n-Ydim2, alpha, A2, lda, X, incX, beta, C2, incY, H2);}
                            );
                        }
                        )
        }
        return Y;
    }


    template<class Field, class AlgoT, class FieldTrait, class Cut>
    typename Field::Element_ptr
    fgemv(const Field& F,
           const FFLAS_TRANSPOSE ta,
           const size_t m,
           const size_t n,
           const typename Field::Element alpha,
           const typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::ConstElement_ptr X, const size_t incX,
           const typename Field::Element beta,
           typename Field::Element_ptr Y, const size_t incY,
           MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, Cut> > & H){
        size_t Xdim, Ydim;
        if (ta == FflasNoTrans) { Xdim = n; Ydim = m;}
        else {Xdim = m; Ydim = n;}
        SYNCH_GROUP(
                    FORBLOCK1D(iter,Ydim,H.parseq,
                               TASK(MODE( READ (A[iter.begin()*((ta==FflasNoTrans)?lda:1)],X) CONSTREFERENCE(F) READWRITE (Y[iter.begin()*incY]) ),
                                    {
                                    if (ta==FflasNoTrans)
                                    fgemv( F, ta, (iter.end()-iter.begin()), Xdim, alpha, A + iter.begin()*lda, lda, X, incX, beta, Y+iter.begin()*incY, incY);
                                    else
                                    fgemv( F, ta, Xdim, (iter.end()-iter.begin()), alpha, A + iter.begin(), lda, X, incX, beta, Y+iter.begin()*incY, incY);
                                    }
                                   )
                              );
                   );
        return Y;
    }


} // FFLAS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
