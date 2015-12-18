/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfgemm.inl
 * Copyright (C) 2015 Jean Guillaume Dumas Clement Pernet Ziad Sultan
 *<ziad.sultan@imag.fr>
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

#include "fflas-ffpack/paladin/parallel.h"

namespace FFLAS
{

	template<class Field>
	void pfzero(const Field& F, 		   
				size_t m, size_t n, 
				typename Field::Element_ptr C, 
				size_t BS=0)
	{
		BS=std::max(BS, (size_t)Protected::WinogradThreshold(F) );
		SYNCH_GROUP(
			for(size_t p=0; p<m; p+=BS) ///row
			for(size_t pp=0; pp<n; pp+=BS) //column
		{
			size_t M=BS, MM=BS;
			if(!(p+BS<m))
				M=m-p;
			if(!(pp+BS<n))
				MM=n-pp;
			TASK(MODE(CONSTREFERENCE(F)),
			{
				for(size_t j=0; j<M; j++)
					for(size_t jj=0; jj<MM; jj++)
						F.assign(C[(p+j)*n+pp+jj],F.zero);
			});
		}
			);
	}

	template<class Field>
	void pfrand(const Field &F,
				size_t m, size_t n, 
				typename Field::Element_ptr C, 
				size_t BS=0)
	{
		BS=std::max(BS, (size_t)Protected::WinogradThreshold(F) );
		typename Field::RandIter G(F); 
		SYNCH_GROUP(
			for(size_t p=0; p<m; p+=BS) ///row
			for(size_t pp=0; pp<n; pp+=BS) //column
		{
			size_t M=BS, MM=BS;
			if(!(p+BS<m))
				M=m-p;
			if(!(pp+BS<n))
				MM=n-pp;
			TASK(MODE(CONSTREFERENCE(G)),
			{
				for(size_t j=0; j<M; j++)
					for(size_t jj=0; jj<MM; jj++)
						G.random (*(C+(p+j)*n+pp+jj));
			});
		}
			);	
	}



} // FFLAS
