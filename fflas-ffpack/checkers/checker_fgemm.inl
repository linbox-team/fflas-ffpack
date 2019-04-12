/* checkers/checker_fgemm.inl
 * Copyright (C) 2016 Ashley Lesdalons
 *
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
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

#ifndef __FFLASFFPACK_checker_fgemm_INL
#define __FFLASFFPACK_checker_fgemm_INL

namespace FFLAS {

    template <class Field>
    class CheckerImplem_fgemm {

        const Field& F;
        const size_t m,n,k,ldc;
        typename Field::Element_ptr v,w1;

    public:
        CheckerImplem_fgemm(const Field &F_,
                            const size_t m_, const size_t n_, const size_t k_,
                            const typename Field::Element beta,
                            typename Field::Element_ptr C, const size_t ldc_)
        : F(F_), m(m_), n(n_), k(k_), ldc(ldc_), v(FFLAS::fflas_new(F_,n)),w1(FFLAS::fflas_new(F_,m))
        {
            typename Field::RandIter G(F);
            init(G,beta,C);
        }

        CheckerImplem_fgemm(typename Field::RandIter &G,
                            const size_t m_, const size_t n_, const size_t k_,
                            const typename Field::Element beta,
                            typename Field::Element_ptr C, const size_t ldc_)
        : F(G.ring()), m(m_), n(n_), k(k_), ldc(ldc_), v(FFLAS::fflas_new(F,n)),w1(FFLAS::fflas_new(F,m))
        {
            init(G,beta,C);
        }

        ~CheckerImplem_fgemm() {
            FFLAS::fflas_delete(v,w1);
        }

        inline bool check(const FFLAS::FFLAS_TRANSPOSE ta,
                          const FFLAS::FFLAS_TRANSPOSE tb,
                          const typename Field::Element alpha,
                          typename Field::ConstElement_ptr A, const size_t lda,
                          typename Field::ConstElement_ptr B, const size_t ldb,
                          typename Field::ConstElement_ptr C)
        {
            // w1 <- C.v - w1
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, C, ldc, v, 1, F.mOne, w1, 1);

            size_t Arows, Acols, Brows, Bcols;
            if (tb == FFLAS::FflasNoTrans){ Brows = k; Bcols = n;} else { Brows = n; Bcols = k;}
            if (ta == FFLAS::FflasNoTrans){ Arows = m; Acols = k;} else { Arows = k; Acols = m;}

            // w2 <- B.v
            typename Field::Element_ptr w2 = FFLAS::fflas_new(F,k);
            FFLAS::finit(F,k,w2,1);
            FFLAS::fgemv(F, tb, Brows, Bcols, F.one, B, ldb, v, 1, F.zero, w2, 1);

            // w1 <- alpha.A.w2 - w1
            FFLAS::fgemv(F, ta, Arows, Acols, alpha, A, lda, w2, 1, F.mOne, w1, 1);

            FFLAS::fflas_delete(w2);

            // is w1 == O ?
            bool pass = FFLAS::fiszero(F, m, w1, 1);
            if (!pass) throw FailureFgemmCheck();
            return pass;
        }

    private:
        inline void init(typename Field::RandIter &G, const typename Field::Element beta, typename Field::Element_ptr C) {
            FFLAS::finit(F,n,v,1);
            FFLAS::finit(F,m,w1,1);
            FFLAS::frand(F,G,n,v,1);

            // w1 <- beta.C.v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, beta, C, ldc, v, 1, F.zero, w1, 1);
        }

    };
}
#endif // __FFLASFFPACK_checker_fgemm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
