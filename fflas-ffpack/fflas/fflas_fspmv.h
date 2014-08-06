/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
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

/** @file fflas/fflas_fspmv.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_fspmv_H
#define __FFLASFFPACK_fflas_fflas_fspmv_H


#ifndef index_t
#define index_t size_t
#endif

#include "fflas-ffpack/config.h"
#include "fflas-ffpack/config-blas.h"

#ifdef __FFLASFFPACK_HAVE_MKL
#ifndef _MKL_H_ // temporary
#error "MKL (mkl.h) not present, while you have MKL enabled"
#endif
#undef index_t
#define index_t MKL_INT
#endif

namespace FFLAS { /*  DNS */

	template<class Element>
	struct VECT ;

	template<class Element>
	struct DNS ;

}

namespace FFLAS { /*  COO */

	template<class Element>
	struct COO ;

	template<class Element>
	struct COO_sub ;

	template<class Element>
	struct COO_ZO ;


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );

	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      const COO_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     );

} // FFLAS

namespace FFLAS { /*  CSR */

	template<class Element>
	struct CSR ;

	template<class Element>
	struct CSR_sub ;

	template<class Element>
	struct CSR_ZO ;


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );


	// y = A.x + b y
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );

	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     );


}

namespace FFLAS { /*  CSC */

} // FFLAS

namespace FFLAS { /*  ELL */

	template<class Element>
	struct ELL ;

	template<class Element>
	struct ELLR ;
} // FFLAS

namespace FFLAS { /* SELL */

	template<class Element>
	struct SELL;

	template<class Element>
	struct SELL_sub;

	template<class Element>
	struct SELL_ZO;

	template<class Field>
	void sp_spmv(const Field & F,
	 			 const SELL<typename Field::Element> & A,
	  			 const VECT<typename Field::Element> & x,
	  			 const typename Field::Element b,
	  			 VECT<typename Field::Element> & y);

	template<class Field>
	void sp_spmv(const Field & F,
	 			 const SELL_sub<typename Field::Element> & A,
	  			 const VECT<typename Field::Element> & x,
	  			 const typename Field::Element b,
	  			 VECT<typename Field::Element> & y);

	template<class Field>
	void sp_spmv(const Field & F,
	 			 const SELL_ZO<typename Field::Element> & A,
	  			 const VECT<typename Field::Element> & x,
	  			 const typename Field::Element b,
	  			 VECT<typename Field::Element> & y);

}

namespace FFLAS { /*  BCSR */

} // FFLAS

namespace FFLAS { /*  DIA */

} // FFLAS

namespace FFLAS { /*  SKY */

} // FFLAS

namespace FFLAS { /*  JAG */

} // FFLAS


#include "fflas-ffpack/fflas/fflas_fspmv.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/fflas_fspmv_coo.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/fflas_fspmv_csr.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/fflas_fspmv_ell.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/fflas_fspmv_ellr.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/fflas_fspmv_sell.inl"

#endif // __FFLASFFPACK_fflas_fflas_fspmv_H
