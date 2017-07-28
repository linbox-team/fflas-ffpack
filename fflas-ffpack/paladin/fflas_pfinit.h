/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_pfinit.inl
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

#ifndef __FFLASFFPACK_pfinit_H
#define __FFLASFFPACK_pfinit_H

#include "fflas-ffpack/paladin/parallel.h"

namespace FFLAS
{

	template<class Field>
	void pfzero(const Field& F, 		   
				size_t m, size_t n, 
				typename Field::Element_ptr C, 
				size_t BS=0)
	{
		using FFLAS::CuttingStrategy::Block;
		using FFLAS::StrategyParameter::Grain;

		BS=std::max(BS, (size_t)Protected::WinogradThreshold(F) );

		SYNCH_GROUP(
            FORBLOCK2D(iter, m, n, SPLITTER(BS, Block, Grain),
                       TASK(MODE(CONSTREFERENCE(F)),
                       {
                           fzero(F, 
                                 iter.iend()-iter.ibegin(),
                                 iter.jend()-iter.jbegin(),
                                 C+iter.ibegin()*n+iter.jbegin(),
                                 n);
                       }
                            );
                       );
            );
	}

	template<class Field, class RandIter>
	void pfrand(const Field& F, 
                RandIter& G,
				size_t m, size_t n, 
				typename Field::Element_ptr C, 
				size_t BS=0)
	{
		using FFLAS::CuttingStrategy::Block;
		using FFLAS::StrategyParameter::Grain;

		BS=std::max(BS, (size_t)Protected::WinogradThreshold(F) );
		SYNCH_GROUP(
            FORBLOCK2D(iter, m, n, SPLITTER(BS, Block, Grain),
                       TASK(MODE(CONSTREFERENCE(F,G)),
                       {
                           frand(F, G, 
                                 iter.iend()-iter.ibegin(),
                                 iter.jend()-iter.jbegin(),
                                 C+iter.ibegin()*n+iter.jbegin(),
                                 n);
                       }
                            );
                       );
			);	
	}


		// d <- d + <x,y>
	template<class Field>
	typename Field::Element& pfdot(const Field& F, const size_t N, 
								   typename Field::ConstElement_ptr x, const size_t incx,
								   typename Field::ConstElement_ptr y, const size_t incy,
								   typename Field::Element& d, size_t rs = NUM_THREADS)
	{
		typename Field::Element * z = fflas_new<typename Field::Element>(rs);
		for(size_t i=0;i<rs;++i) F.init(z[i]);

		SYNCH_GROUP({
			FORBLOCK1D(iter, N, SPLITTER(rs), {
				const size_t i( iter.begin()/(iter.end()-iter.begin()) );					
				TASK(MODE(CONSTREFERENCE(F) READ(x,y) WRITE(z[i])), {
// 						std::cerr << NUM_THREADS << ", T(" << omp_get_thread_num()<< "):" << iter.end()-iter.begin() << std::endl;
// 						for(auto itee=iter.begin(); itee!=iter.end(); ++itee)
// 							F.write(F.write(std::cerr, *(x+itee)) << '.', *(y+itee)) << std::endl;
					F.assign(z[i], fdot(F, 
										iter.end()-iter.begin(),
										x+iter.begin()*incx, incx,
										y+iter.begin()*incy, incy));
// 						F.write(std::cerr << "pd:", z) << std::endl;
				});
			});
		});

		for(size_t i=0;i<rs;++i) F.addin(d,z[i]);
		return d;
	}

 
	template<class Field>
	typename Field::Element pfdot(const Field& F, const size_t N, 
								  typename Field::ConstElement_ptr x, const size_t incx,
								  typename Field::ConstElement_ptr y, const size_t incy,
								  size_t rs = NUM_THREADS)
	{
		typename Field::Element d; F.init(d); F.assign(d,F.zero);
		
		PAR_BLOCK {
			pfdot(F, N, x, incx, y, incy, d, rs);
		}
		
		return d;
	}



} // FFLAS

#endif
