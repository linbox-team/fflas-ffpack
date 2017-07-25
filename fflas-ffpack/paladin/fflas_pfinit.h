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

	template<class Field>
	typename Field::Element pfdot(const Field& F, const size_t N, 
								   typename Field::ConstElement_ptr x, const size_t incx,
								   typename Field::ConstElement_ptr y, const size_t incy)
	{
		typename Field::Element d; F.init(d); F.assign(d,F.zero);
		
		PAR_BLOCK {
			
			SYNCH_GROUP(
				FORBLOCK1D(iter, N, SPLITTER(),
						   TASK(MODE(CONSTREFERENCE(F) READ(x,y) READWRITE(d)),
						   {
							   std::cerr << NUM_THREADS << ", T(" << omp_get_thread_num()<< "):" << iter.end()-iter.begin() << std::endl;
//                       for(auto itee=iter.begin(); itee!=iter.end(); ++itee)
//                           F.write(F.write(std::cerr, *(x+itee)) << '.', *(y+itee)) << std::endl;
							   
							   typename Field::Element z; F.init(z); F.assign(z,F.zero); 
							   F.assign(z, fdot(F, 
												iter.end()-iter.begin(),
												x+iter.begin(), incx,
												y+iter.begin(), incy));	
#pragma omp critical
							   F.addin(d,z);
							   
						   }
								);
						   );
				WAIT;
				);
//         typename Field::Element d; F.init(d); F.assign(d,F.zero);
//         for(size_t i=0; i<rs; ++i) {
//             F.write(std::cerr << "z[" << i << "]:", z[i]) << std::endl;
//             F.addin(d,z[i]);
//         }
//         return d;
			
		}
        return d;
	}



} // FFLAS

#endif
