/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/Checker_charpoly.inl
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

#ifndef __FFLASFFPACK_checker_charpoly_INL
#define __FFLASFFPACK_checker_charpoly_INL

#ifdef ENABLE_CHECKER_charpoly
#include "ffpack/ffpack.h"

namespace FFPACK {
    template <class Field, class Polynomial> 
    class Checker_charpoly {

        const Field& F;
        const size_t n;
        typename Field::Element lambda, det;
        bool pass;

    public:
        Checker_charpoly(const Field& F_, const size_t n_, typename Field::Element_ptr A) 
                : F(F_), n(n_)
            {
                typename Field::RandIter G(F);
                init(G,A);
            }

        Checker_charpoly(typename Field::RandIter &G, const size_t n_, typename Field::Element_ptr A)
                : F(G.ring()), n(n_)
            {
                init(G,A);
            }

        ~Checker_charpoly() {
        }

        inline bool check(Polynomial &g) {
            typename Field::Element h = F.zero,
                t = F.one,
                u;
            for (size_t i=0; i < g.size(); ++i) {
                F.mul(u,g[i],t);
                F.add(h,h,u);
                F.mul(t,t,lambda);
            }

                // is h == det ?
            pass = pass && F.areEqual(h,det);
            if (!pass) throw FailureCharpolyCheck();

            return pass;
        }

    private:
        inline void init(typename Field::RandIter &G, typename Field::Element_ptr A) {
                // random lambda
            G.random(lambda);

            typename Field::Element_ptr v = FFLAS::fflas_new(F,n,1),
                w = FFLAS::fflas_new(F,n,1),
                Ac = FFLAS::fflas_new(F,n,n);
            FFLAS::frand(F,G,n,v,1);

                // w <- -A.v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.mOne, A, n, v, 1, F.zero, w, 1);

            if (!F.isZero(lambda)) {
                    // w <- lambda.v + w
                FFLAS::faxpy(F, n, lambda, v, 1, w, 1);
            }

                // Ac <- A - lambda.I
            for (size_t i=0; i<n; ++i)
                for (size_t j=0; j<n; ++j) {
                    if (i==j) F.sub(*(Ac+i*n+i),*(A+i*n+i),lambda);
                    else F.assign(*(Ac+i*n+j),*(A+i*n+j));
                }

                // w <- Ac.v + w
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, Ac, n, v, 1, F.one, w, 1);

                // is w == 0 ?
            pass = FFLAS::fiszero(F,n,1,w,1);
            FFLAS::fflas_delete(v,w);
            if (!pass) throw FailureCharpolyCheck();

                // P,Ac,Q <- PLUQ(Ac)
            size_t *P = FFLAS::fflas_new<size_t>(n);
            size_t *Q = FFLAS::fflas_new<size_t>(n);
            FFPACK::PLUQ(F, FFLAS::FflasNonUnit, n, n, Ac, n, P, Q);

                // compute the determinant of A
            F.init(det,*Ac);
            for (size_t i=1; i<n; ++i)
                F.mul(det,det,*(Ac+i*n+i));
            if (n%2 == 1) F.neg(det,det);

                // count the number of permutations
            int t = 0;
            for (size_t i=0; i<n; ++i) {
                if (P[i] != i) t++;
                if (Q[i] != i) t++;
            }
            if (t%2 == 1) F.neg(det,det);

            FFLAS::fflas_delete(Ac);
        }
    };
    
}

#endif // ENABLE_CHECKER_charpoly

#endif // __FFLASFFPACK_checker_charpoly_INL
