/* checkers/checker_det.inl
 * Copyright (C) 2016 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __FFLASFFPACK_checker_det_INL
#define __FFLASFFPACK_checker_det_INL

#include "fflas-ffpack/ffpack/ffpack.h"

#ifdef TIME_CHECKER_Det
#include <givaro/givtimer.h>
#endif

namespace FFPACK {

    // Checker for the determinant following [Dumas Lucas Pernet'17].
    // First implementation, meant to evolve into an RPM checker and/or
    // and interractive Det checker.
    template <class Field>
    class CheckerImplem_Det {

        const Field& F;
        typename Field::Element_ptr u,v,w;
        typename Field::Element du,dv;
        const size_t n;
#ifdef TIME_CHECKER_Det
        mutable Givaro::Timer _time;
#endif

    public:
        CheckerImplem_Det(const Field& F_, size_t n_,
                          typename Field::ConstElement_ptr A, size_t lda)
        : F(F_),
        u(FFLAS::fflas_new(F_,n_,2)),
        v(u+1),
        w(FFLAS::fflas_new(F_,n_)),
        n(n_)
        {
            typename Field::RandIter G(F);
            init(G,A,lda);
        }

        CheckerImplem_Det(typename Field::RandIter &G, size_t n_,
                          typename Field::ConstElement_ptr A, size_t lda)
        : F(G.ring()),
        u(FFLAS::fflas_new(F,n_,2)),
        v(u+1),
        w(FFLAS::fflas_new(F,n_)),
        n(n_)
        {
            init(G,A,lda);
        }

        ~CheckerImplem_Det() {
            FFLAS::fflas_delete(u,w);
        }

        /** check if the Det factorization is correct.
         *  Needs matrix in LU form
         * @param LU, storage for L and U
         * @param det
         * @param P
         * @param Q
         */
        inline bool check(const typename Field::Element& det,
                          typename Field::ConstElement_ptr LU, size_t lda,
                          size_t *P, size_t *Q) const {
#ifdef TIME_CHECKER_Det
            Givaro::Timer checktime, overhead; checktime.start();
#endif
            // u <-- Q.u, v <-- Q.v
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 2, 0, n, u, 2, Q);

            // w <-- w P
            FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasTrans, 1, 0, n, w, 1, P);

#ifdef TIME_CHECKER_Det
            checktime.stop(); _time += checktime;
            overhead.start();
#endif
            // u <-- U.u, v <-- U.v
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, n, 2, F.one, LU, lda, u, 2);

            // w <-- w.L
            // Warning: should be ftrmv
            FFLAS::ftrmm(F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, 1, n, F.one, LU, lda, w, n);

#ifdef TIME_CHECKER_Det
            overhead.stop();
            checktime.start();
#endif
            typename Field::Element zu, zv;
            zu = FFLAS::fdot(F, n, w, 1, u, 2);
            zv = FFLAS::fdot(F, n, w, 1, v, 2);

            bool pass = F.areEqual(zu,du) && F.areEqual(zv,dv);
            if (!pass) throw FailureDetCheck();

            // Check det
            typename Field::Element dd;
            F.init (dd);
            F.assign (dd,F.one);

            for (typename Field::ConstElement_ptr Ai(LU); Ai < LU+ n*lda+n; Ai+=lda+1 )
                F.mulin( dd, *Ai );

            int count=0;
            for (size_t i=0;i<n;++i)
                if (P[i] != i) ++count;
            for (size_t i=0;i<n;++i)
                if (Q[i] != i) ++count;
            if ((count&1) == 1)
                F.negin(dd);

            pass &= F.areEqual(dd,det);
            if (!pass) throw FailureDetCheck();


#ifdef TIME_CHECKER_Det
            checktime.stop(); _time += checktime;
            std::cout << "Det OVERH: " << overhead << std::endl;
            std::cout << "Det CHECK: " << _time << std::endl;
#endif
            return pass;
        }

    private:
        inline void init(typename Field::RandIter &G,
                         typename Field::ConstElement_ptr A, size_t lda) {
#ifdef TIME_CHECKER_Det
            Givaro::Timer inittime; inittime.start();
#endif
            FFLAS::finit(F,n*2,u,1);
            FFLAS::finit(F,n,w,1);
            FFLAS::frand(F,G,n,2,u,2);
            FFLAS::frand(F,G,n,w,1);

            typename Field::Element_ptr t(FFLAS::fflas_new(F,n));
            FFLAS::finit(F,n,t,1);

            // t <-- w . A
            FFLAS::fgemv(F, FFLAS::FflasTrans, n, n, F.one, A, lda, w, 1, F.zero, t, 1);
            du = FFLAS::fdot(F, n, u, 2, t, 1);
            dv = FFLAS::fdot(F, n, v, 2, t, 1);

            FFLAS::fflas_delete(t);
#ifdef TIME_CHECKER_Det
            inittime.stop(); _time += inittime;
#endif
        }
    };
}
#endif // __FFLASFFPACK_checker_det_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
