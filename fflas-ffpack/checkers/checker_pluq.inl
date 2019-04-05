/* checkers/checker_pluq.inl
 * Copyright (C) 2016 Jean-Guillaume Dumas
 *
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_checker_pluq_INL
#define __FFLASFFPACK_checker_pluq_INL

#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/fflas_io.h"

#ifdef TIME_CHECKER_PLUQ
#include <givaro/givtimer.h>
#endif

namespace FFPACK {
    template <class Field>
    class CheckerImplem_PLUQ {

        const Field& F;
        typename Field::Element_ptr v,w;
        const size_t m,n;
#ifdef TIME_CHECKER_PLUQ
        mutable Givaro::Timer _time;
#endif

    public:
        CheckerImplem_PLUQ(const Field& F_, size_t m_, size_t n_,
                           typename Field::ConstElement_ptr A, size_t lda)
        : F(F_),
        v(FFLAS::fflas_new(F_,n_)),
        w(FFLAS::fflas_new(F_,m_)),
        m(m_), n(n_)
        {
            typename Field::RandIter G(F);
            init(G,A,lda);
        }

        CheckerImplem_PLUQ(typename Field::RandIter &G, size_t m_, size_t n_,
                           typename Field::ConstElement_ptr A, size_t lda)
        : F(G.ring()),
        v(FFLAS::fflas_new(F,n_)),
        w(FFLAS::fflas_new(F,m_)),
        m(m_), n(n_)
        {
            init(G,A,lda);
        }

        ~CheckerImplem_PLUQ() {
            FFLAS::fflas_delete(v,w);
        }

        /** check if the PLUQ factorization is correct.
         *  Returns true if w - P(L(U(Q.v))) == 0
         * @param A
         * @param r
         * @param P
         * @param Q
         */
        inline bool check(typename Field::ConstElement_ptr A, size_t lda,
                          const FFLAS::FFLAS_DIAG Diag,
                          size_t r, size_t *P, size_t *Q) const {
#ifdef TIME_CHECKER_PLUQ
            Givaro::Timer checktime; checktime.start();
#endif
            // _w = [w1|w2]
            typename Field::Element_ptr _w = FFLAS::fflas_new(F,std::max(m,n));

            // _w <- v
            // WARNING check fassign on vectors
            FFLAS::fassign(F, n, 1, v, 1, _w, 1);

            // _w <-- Q._w
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, _w, 1, Q);


            // w1 <- U1.w1
            // WARNING: could be ftrmv
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, r, 1, F.one, A, lda, _w, 1);

            // w1 <- U2.w2 + w1
            if (r < n)
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r, 1, n-r, F.one, A+r, lda, _w+r, 1, F.one, _w, 1);

            // w2 <- L2.w1
            if (r < m)
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m-r, 1, r, F.one, A+r*lda, lda, _w, 1, F.zero, _w+r, 1);

            const FFLAS::FFLAS_DIAG oppDiag = (Diag == FFLAS::FflasNonUnit) ? FFLAS::FflasUnit : FFLAS::FflasNonUnit;
            // w1 <- L1.w1
            // WARNING: should be ftrmv
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, oppDiag, r, 1, F.one, A, lda, _w, 1);

            // _w <- P._w
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans, 1, 0, m, _w, 1, P);

            // is _w == w ?
            bool pass = FFLAS::fequal(F,m,w,1,_w,1);
            FFLAS::fflas_delete(_w);

            if (!pass) throw FailurePLUQCheck();

#ifdef TIME_CHECKER_PLUQ
            checktime.stop(); _time += checktime;
            std::cerr << "PLUQ CHECK: " << _time << std::endl;
#endif
            return pass;
        }

    private:
        inline void init(typename Field::RandIter &G,
                         typename Field::ConstElement_ptr A, size_t lda) {
#ifdef TIME_CHECKER_PLUQ
            Givaro::Timer inittime; inittime.start();
#endif
            FFLAS::finit(F,n,v,1);
            FFLAS::finit(F,m,w,1);
            FFLAS::frand(F,G,n,v,1);

            // w <-- A.v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, lda, v, 1, F.zero, w, 1);
#ifdef TIME_CHECKER_PLUQ
            inittime.stop(); _time += inittime;
#endif
        }
    };
}
#endif // __FFLASFFPACK_checker_pluq_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
