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

/** @file fflas/fflas_fspmv.inl
*/

#ifndef __FFLASFFPACK_fflas_fflas_fspmv_INL
#define __FFLASFFPACK_fflas_fflas_fspmv_INL

#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/fflas/fflas_bounds.inl"

namespace FFLAS { /*  DNS */

	template<class Element>
	struct VECT {
		size_t m ;
		size_t inc ;
		Element * dat ;
	};

	template<class Element>
	struct DNS {

		size_t n ;
		size_t ld ;
		Element * dat ;
	};

} // FFLAS

namespace FFLAS { /* HYB */
#if 0
	template<class Element>
	struct SPADD {
		size_t ncsr;
		CSR<Element> * csr;
		size_t ncoo;
		COO<Element> * coo;
		size_t ndns;
		DNS<Element> * dns;
		size_t nell;
		ELL<Element> * ell;
		size_t nellr ;
		ELLR<Element> * ellr ;
		size_t ndia ;
		DIA<Element> * dia;

		SPADD() :
			ncsr(0)  ,csr(NULL)
			,ncoo(0) ,coo(NULL)
			,ndns(0) ,dns(NULL)
			,ndia(0) ,dia(NULL)
			,nell(0) ,ell(NULL)
			,nellr(0),ellr(NULL)
		{}
	};
#endif

} // FFLAS


namespace FFLAS { /*  BCSR */

} // FFLAS

namespace FFLAS { /*  DIA */

} // FFLAS

namespace FFLAS { /*  SKY */

} // FFLAS

namespace FFLAS { /*  JAG */

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_fspmv_INL
