/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
 *              Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_fspmv_ell.inl
 * NO DOC
*/

#ifndef __FFLASFFPACK_fflas_fflas_spmv_ell_INL
#define __FFLASFFPACK_fflas_fflas_spmv_ell_INL

namespace FFLAS { /*  ELL */

	template<class Element>
	struct ELL {
		size_t m ;
		size_t n ;
		size_t  ld ;
		index_t  * col ;
		Element * dat ;
	};

	template<class Element>
	struct ELL_sub : public ELL<Element> {
	}

	namespace details {

	// y = A x + b y ; (generic)
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const index_t * col,
			      const size_t ld,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      const typename Field::Element & b,
			      typename Field::Element * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if (! F.isOne(b)) {
					if (F.isZero(b)) {
						F.assign(y[i],F.zero);
					}
					else if (F.isMone(b)) {
						F.negin(y[i]);
					}
					else {
						F.mulin(y[i],b);
					}
				}
				// XXX can be delayed
				for (index_t j = 0 ; j < ld ; ++j) {
					if (F.isZero(dat[i*ld+j]))
						break;
					F.axpyin(y[i],dat[i*ld+j],x[col[i*ld+j]]);
				}
			}
		}

		// double
		template<>
		void sp_fgemv(
			      const DoubleDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const index_t * col,
			      const size_t ld,
			      const double*  dat,
			      const double* x ,
			      const double& b,
			      double * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if ( b != 1) {
					if ( b == 0.) {
						y[i] = 0;
					}
					else if ( b == -1 ) {
						y[i]= -y[i];
					}
					else {
						y[i] = y[i] * b;
					}
				}
				for (index_t j = 0 ; j < ld ; ++j) {
					if (dat[i*ld+j] == 0)
						break;
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}

		// float
		template<>
		void sp_fgemv(
			      const DoubleDomain& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const index_t * col,
			      const size_t ld,
			      const float*  dat,
			      const float* x ,
			      const float& b,
			      float * y
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				if ( b != 1) {
					if ( b == 0.) {
						y[i] = 0;
					}
					else if ( b == -1 ) {
						y[i]= -y[i];
					}
					else {
						y[i] = y[i] * b;
					}
				}
				for (index_t j = 0 ; j < ld ; ++j) {
					if (dat[i*ld+j] == 0)
						break;
					y[i] += dat[i*ld+j]*x[col[i*ld+j]];
				}
			}
		}


		// delayed by kmax
		template<class Field>
		void sp_fgemv(
			      const Field& F,
			      // const FFLAS_TRANSPOSE tA,
			      const size_t m,
			      const size_t n,
			      const index_t * col,
			      const size_t ld,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      // const typename Field::Element & b,
			      typename Field::Element * y,
			      const index_t & kmax
			     )
		{
			for (size_t i = 0 ; i < m ; ++i) {
				// y[i] = 0;
				index_t j = st[i];
				index_t j_loc = j;
				index_t j_end = st[i+1];
				index_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j) {
						y[i] += dat[j] * x[col[j]];
					}
					F.init(y[i],y[i]);
				}
				for ( ; j < j_end ; ++j) {
					y[i] += dat[j] * x[col[j]];
				}
				F.init(y[i],y[i]);
			}
		}
	} // details

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_spmv_ell_INL

