/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/Checker_invert.inl
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

#ifndef __FFLASFFPACK_checker_invert_INL
#define __FFLASFFPACK_checker_invert_INL

#ifdef ENABLE_CHECKER_invert

namespace FFPACK {

    template <class Field> 
    class Checker_invert {

        const Field& F;
        typename Field::Element_ptr v,w;
        const size_t m,lda;

    public:
        Checker_invert(const Field& F_, const size_t m_, typename Field::Element_ptr A, const size_t lda_) 
                : F(F_), v(FFLAS::fflas_new(F_,m_,1)), w(FFLAS::fflas_new(F_,m_,1)), m(m_), lda(lda_)
            {
                typename Field::RandIter G(F);
                init(G,m,A,lda);
            }

        Checker_invert(typename Field::RandIter &G, const size_t m_, typename Field::Element_ptr A, const size_t lda_) 
                : F(G.ring()), v(FFLAS::fflas_new(F,m_,1)), w(FFLAS::fflas_new(F,m_,1)), m(m_), lda(lda_)
            {
                init(G,m,A,lda);
            }

        ~Checker_invert() {
            FFLAS::fflas_delete(v,w);
        }

        inline bool check(typename Field::Element_ptr A, int nullity) {
                // v <- A.w - v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, w, 1, F.mOne, v, 1);

            bool pass = FFLAS::fiszero(F,m,1,v,1) || (nullity != 0);
            if (!pass) throw FailureInvertCheck();
            return pass;
        }

    private:
        void init(typename Field::RandIter &G, const size_t m_, typename Field::Element_ptr A, const size_t lda_) {
            FFLAS::frand(F,G,m,v,1);

// write_field(F,std::cerr<<"init A : ",A,m,m,lda,true)<<std::endl;
// write_field(F,std::cerr<<"init v : ",v,m,1,1,true)<<std::endl;  
    	
                // w <- A.v
            FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, m, F.one, A, lda, v, 1, F.zero, w, 1);
// write_field(F,std::cerr<<"init w : ",w,m,1,1,true)<<std::endl;
        }
  
    };
}
#endif // ENABLE_CHECKER_invert

#endif // __FFLASFFPACK_checker_invert_INL
