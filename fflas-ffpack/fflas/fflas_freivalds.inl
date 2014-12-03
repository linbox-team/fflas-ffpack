/* fflas/fflas_freivalds.inl
 * Copyright (C) 2014 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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
#ifndef __FFLASFFPACK_freivalds_INL
#define __FFLASFFPACK_freivalds_INL

// #include "fflas-ffpack/utils/Matio.h"

namespace FFLAS{ 

	/** @brief  freivalds: <b>F</b>reivalds <b>GE</b>neral <b>M</b>atrix <b>M</b>ultiply <b>R</b>andom <b>C</b>heck.
	 *
	 * Randomly Checks \f$C = \alpha \mathrm{op}(A) \times \mathrm{op}(B)\f$
	 * \param F field.
	 * \param ta if \c ta==FflasTrans then \f$\mathrm{op}(A)=A^t\f$, else \f$\mathrm{op}(A)=A\f$,
	 * \param tb same for matrix \p B
	 * \param m see \p A
	 * \param n see \p B
	 * \param k see \p A
	 * \param alpha scalar
	 * \param A \f$\mathrm{op}(A)\f$ is \f$m \times k\f$
	 * \param B \f$\mathrm{op}(B)\f$ is \f$k \times n\f$
	 * \param C \f$C\f$ is \f$m \times n\f$
	 * \param lda leading dimension of \p A
	 * \param ldb leading dimension of \p B
	 * \param ldc leading dimension of \p C
	 */
    template<class Field> inline  bool
    freivalds (const Field& F,
               const FFLAS_TRANSPOSE ta,
               const FFLAS_TRANSPOSE tb,
               const size_t m, const size_t n, const size_t k,
               const typename Field::Element alpha,
               typename Field::ConstElement_ptr A, const size_t lda,
               typename Field::ConstElement_ptr B, const size_t ldb,
               typename Field::ConstElement_ptr C, const size_t ldc) {
    
        typename Field::Element_ptr v, y, x;

        v = FFLAS::fflas_new<typename Field::Element>(n);
        y = FFLAS::fflas_new<typename Field::Element>(k);
        x = FFLAS::fflas_new<typename Field::Element>(m);

        typename Field::RandIter G(F);
        for(size_t j=0; j<(size_t)n; ++j)
            G.random(*(v+j));
    
//         F.write(std::cerr<< "alpha:", alpha) << std::endl;
//         write_field(F,std::cerr<<"A:",A,m,k,lda,true) << std::endl;
//         write_field(F,std::cerr<<"B:",B,k,n,ldb,true) << std::endl;
//         write_field(F,std::cerr<<"C:",C,m,n,ldc,true) << std::endl;

            // y <-- 1.\mathrm{op}(B).v
        FFLAS::fgemv(F, tb, k,n, 
                     F.one, B, n, v, 1, F.zero, y, 1);
            // x <-- alpha.\mathrm{op}(A).y
            // x <-- alpha.\mathrm{op}(A).\mathrm{op}(B).v
        FFLAS::fgemv(F, ta, m,k, 
                     alpha, A, k, y, 1, F.zero, x, 1);
 
            // x <-- -C.v+x =?= 0
        FFLAS::fgemv(F, FFLAS::FflasNoTrans,m,n, 
                     F.mOne, C, n, v, 1, F.one, x, 1);

    
        bool pass=true;
        for(size_t j=0; j<(size_t)m; ++j) 
            pass &= F.isZero( *(x+j) );

        return pass;
    }
 
}

    

#endif //  __FFLASFFPACK_freivalds_INL
    
    
