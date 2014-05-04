/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by BB <bbboyer@ncsu.edu>
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

/** @file fflas_fgemm/fgemm_winograd.h
 * @brief Strassen--Winograd matrix multiplication.
 * @warning The domain is supposed to be a field since some divisions are required for efficiency purposes
 * An alternative has to be written for finite rings if necessary
 */

#ifndef __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL
#define __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL

#include "fflas_bounds_winograd.inl"

namespace FFLAS {

	struct Winograd2Helper : public MMParameters {
		int w ;
		size_t kmax ;
		FFLAS_BASE base ;

		Winograd2Helper() :
			w(-1) , kmax(0), base(FflasDouble)
		{}

		Winograd2Helper(int rec) :
			w(rec) , kmax(0), base(FflasDouble)
		{}

		template<class Field>
		void computeParameters(const Field & F
				  , const size_t & m
				  , const size_t & n
				  , const size_t & k
				  , const typename Field::Element &alpha
				  , const typename Field::Element &beta)
		{
			bool winoLevelProvided = (w != (int(-1)));
			// size_t kmax = 0;
			int winolevel = w;
			// FFLAS_BASE base;
			typename Field::Element gamma;
			F.div(gamma,beta,alpha);
			Protected::MatMulParametersWinograd (F, m, n, k, gamma, kmax, base,
						     winolevel, winoLevelProvided);
			w = winolevel;

		}
	};
}

namespace FFLAS { namespace Protected {

} // Protected

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_fgemm_winograd_INL
