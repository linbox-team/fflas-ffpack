/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
    template <class Field> 
    class CheckerImplem_Det {

        const Field& F;
        typename Field::Element_ptr u,v,w;
		typename Field::Element du,dv;
        const size_t n;
#ifdef TIME_CHECKER_Det
        Givaro::Timer _time;
#endif

    public:
        CheckerImplem_Det(const Field& F_, size_t n_, 
                     typename Field::ConstElement_ptr A, size_t lda) 
				: F(F_), 
                  u(FFLAS::fflas_new(F_,n_,1)), 
                  v(FFLAS::fflas_new(F_,n_,1)), 
                  w(FFLAS::fflas_new(F_,n_,1)), 
                  n(n_)
            {
                typename Field::RandIter G(F);
                init(G,A,lda);
            }

        CheckerImplem_Det(typename Field::RandIter &G, size_t n_, 
                     typename Field::ConstElement_ptr A, size_t lda)
				: F(G.ring()), 
                  u(FFLAS::fflas_new(F,n_,1)), 
                  v(FFLAS::fflas_new(F,n_,1)), 
                  w(FFLAS::fflas_new(F,n_,1)), 
                  n(n_)
            {
                init(G,A,lda);
            }

        ~CheckerImplem_Det() {
            FFLAS::fflas_delete(u,v,w);
        }

            /** check if the Det factorization is correct.
             *  Returns true if w - P(L(U(Q.v))) == 0
             * @param A
             * @param r
             * @param P
             * @param Q
             */
        inline bool check(typename Field::ConstElement_ptr A, size_t lda, 
                          const FFLAS::FFLAS_DIAG Diag,
                          size_t *P, size_t *Q) const {
#ifdef TIME_CHECKER_Det
            Givaro::Timer checktime; checktime.start();
#endif
                // u <-- Q.u, v <-- Q.v
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, u, 1, Q);
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);

				// w <-- w.P
            FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, w, 1, P);


				// u <-- U.u, v <-- U.v
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, n, 1, F.one, A, lda, u, 1);
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, n, 1, F.one, A, lda, v, 1);

			typename Field::Element zu, zv;
			zu = FFLAS::fdot(F, n, w, 1, u, 1);
			zv = FFLAS::fdot(F, n, w, 1, v, 1);
			

			bool pass = F.areEqual(zu,du) && F.areEqual(zv,dv);
            if (!pass) throw FailureDetCheck();

#ifdef TIME_CHECKER_Det
            checktime.stop(); _time += checktime;
            std::cerr << "Det CHECK: " << _time << std::endl;
#endif
            return pass;
        }

    private:	
        inline void init(typename Field::RandIter &G, 
                         typename Field::ConstElement_ptr A, size_t lda) {
#ifdef TIME_CHECKER_Det
            Givaro::Timer inittime; inittime.start();
#endif
            FFLAS::frand(F,G,n,u,1);
            FFLAS::frand(F,G,n,v,1);
            FFLAS::frand(F,G,n,w,1);

			typename Field::Element_ptr t(FFLAS::fflas_new(F,n,1));

                // t <-- A.u
			FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, lda, u, 1, F.zero, t, 1);
			du = FFLAS::fdot(F, n, w, 1, t, 1);
			
                // t <-- A.v
			FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, lda, v, 1, F.zero, w, 1);
			dv = FFLAS::fdot(F, n, w, 1, t, 1);
			
            FFLAS::fflas_delete(t);
#ifdef TIME_CHECKER_Det
            inittime.stop(); _time += inittime;
#endif
        }
    };
}
#endif // __FFLASFFPACK_checker_det_INL
