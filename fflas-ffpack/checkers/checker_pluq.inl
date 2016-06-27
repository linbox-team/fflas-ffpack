/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifdef ENABLE_CHECKER_PLUQ
#include "ffpack/ffpack.h"

#ifdef TIME_CHECKER_PLUQ
#include <givaro/givtimer.h>
#endif

namespace FFPACK {
    template <class Field> 
    class Checker_PLUQ {

        const Field& F;
        typename Field::Element_ptr v,w;
        const size_t m,n;
#ifdef TIME_CHECKER_PLUQ
        Givaro::Timer _time;
#endif

    public:
        Checker_PLUQ(const Field& F_, size_t m_, size_t n_, 
                     typename Field::ConstElement_ptr A, size_t lda) 
				: F(F_), 
                  v(FFLAS::fflas_new(F_,n_,1)), 
                  w(FFLAS::fflas_new(F_,m_,1)), 
                  m(m_), n(n_)
            {
                typename Field::RandIter G(F);
                init(G,A,lda);
            }

        Checker_PLUQ(typename Field::RandIter &G, size_t m_, size_t n_, 
                     typename Field::Element_ptr A, size_t lda)
				: F(G.ring()), 
                  v(FFLAS::fflas_new(F,n_,1)), 
                  w(FFLAS::fflas_new(F,m_,1)), 
                  m(m_), n(n_)
            {
                init(G,A,lda);
            }

        ~Checker_PLUQ() {
            FFLAS::fflas_delete(v,w);
        }

            /** check if the PLUQ factorization is correct.
             *  Returns true if w - P(L(U(Q.v))) == 0
             * @param A
             * @param r
             * @param P
             * @param Q
             */
        inline bool check(typename Field::Element_ptr A, size_t lda, 
                          size_t r, size_t *P, size_t *Q) {
#ifdef TIME_CHECKER_PLUQ
            Givaro::Timer checktime; checktime.start();
#endif
				// _w = [w1|w2]
            typename Field::Element_ptr _w = FFLAS::fflas_new(F,m,1); 

                // v <-- Q.v
            FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);

                // w1 <- V1 && w2 <- 0
            FFLAS::fassign(F, r, 1, v, 1, _w, 1);
            FFLAS::fzero(F, m-r, _w+r, 1);

                // w1 <- U1.w1
                // WARNING: should be ftrmv
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, r, 1, F.one, A, lda, _w, 1);
		
                // w1 <- U2.V2 + w1
            if (r < n)
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r, 1, n-r, F.one, A+r, lda, v+r, 1, F.one, _w, 1);

                // w2 <- L2.w1
            if (r < m)
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m-r, 1, r, F.one, A+r*n, lda, _w, 1, F.zero, _w+r, 1); 		

                // w1 <- L1.w1
                // WARNING: should be ftrmv
            FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, 1, F.one, A, lda, _w, 1);

                // _w <- P._w
            FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, _w, 1, P);

				// is _w == w ?
            FFLAS::fsubin(F, m, w, 1, _w, 1);
            bool pass = FFLAS::fiszero(F,m,_w,1);
        
            FFLAS::fflas_delete(_w);

            if (!pass) throw FailurePLUQcheck();

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
            FFLAS::frand(F,G,n,v,1);
    	
                // w <-- A.v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, lda, v, 1, F.zero, w, 1);
#ifdef TIME_CHECKER_PLUQ
            inittime.stop(); _time += inittime;
#endif
        }
    };
}
#endif // ENABLE_CHECKER_PLUQ

#endif // __FFLASFFPACK_checker_pluq_INL
