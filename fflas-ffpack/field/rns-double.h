/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

/*! @file field/rns-double.h
 * @ingroup field
 * @brief  rns structure with double support
 */

#ifndef __FFPACK_rns_double_H
#define __FFPACK_rns_double_H

// Bigger multiple of s lesser or equal than x, s must be a power of two
#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#endif

#include <vector>
#include <givaro/modular-double.h>
#include <givaro/givinteger.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "fflas-ffpack/field/rns-double-elt.h"

namespace FFPACK {

	/* Structure that handles rns representation given a bound and bitsize for prime moduli
	 * support sign representation (i.e. the bound must be twice larger then ||A||)
	 */
	struct rns_double {
		typedef Givaro::Integer integer;
		typedef Givaro::Modular<double> ModField;
		
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basis; // the rns moduli (mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basisMax; // (mi-1)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _negbasis; // (-mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>    _invbasis; // the inverse of rns moduli (1/mi)
		std::vector<ModField> _field_rns; // the associated prime field for each mi
		integer                  _M; // the product of the mi's
		std::vector<integer>         _Mi; // _M/mi
		std::vector<double>         _MMi; // (_Mi)^(-1) mod mi
		std::vector<double>      _crt_in; //  2^(16*j) mod mi
		std::vector<double>     _crt_out; //  (_Mi._MMi) written in base 2^16
		size_t                _size; // the size of the rns basis (number of mi's)
		size_t               _pbits; // the size in bit of the mi's
		size_t                 _ldm; // log[2^16](_M)

		typedef double                        BasisElement;
		typedef rns_double_elt                     Element;
		typedef rns_double_elt_ptr             Element_ptr;
		typedef rns_double_elt_cstptr     ConstElement_ptr;

		rns_double(const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
		:  _M(1), _size(0), _pbits(pbits)
		{
			integer::seeding(seed);
			integer prime;
			integer sum=1;
			while (_M < bound*sum) {
				_basis.resize(_size+1);
				_negbasis.resize(_size+1);
				_basisMax.resize(_size+1);
				do {
					integer::random_exact_2exp(prime, _pbits-1);
					nextprime(prime, prime);
				} while (_M%prime == 0);
				_basis[_size]=prime;
				_basisMax[_size] = prime-1;
				_negbasis[_size] = 0-prime;
				_size++;
				_M*=prime;
				if (rnsmod) sum+=prime;
			}
			precompute_cst();
		}

		template<typename Vect>
		rns_double(const Vect& basis, bool rnsmod=false, long seed=time(NULL))
			:  _basis(basis.begin(),basis.end()), _M(1), _size(basis.size()), _pbits(0)
		{
			for(size_t i=0;i<_size;i++){
				_M*=_basis[i];
				_pbits=std::max(_pbits, integer(_basis[i]).bitsize());
				_basisMax[i] = _basis[i]-1;
				_negbasis[i] = 0-_basis[i];
			}
			std::cout<<"M="<<_M<<std::endl;
			precompute_cst();
		}


		void precompute_cst(){
			_ldm = (_M.bitsize()/16) + ((_M.bitsize()%16)?1:0) ;
			_invbasis.resize(_size);
			_field_rns.resize(_size);
			_Mi.resize(_size);
			_MMi.resize(_size);
			_crt_in.resize(_size*_ldm);
			_crt_out.resize(_size*_ldm);
			const unsigned int MASK=0xFFFF;
			for (size_t i=0;i<_size;i++){
				_invbasis[i]  = 1./_basis[i];
				_field_rns[i] = ModField(_basis[i]);
				_Mi[i]        = _M/(unsigned long)_basis[i];
				_field_rns[i].init(_MMi[i], _Mi[i] % (double)_basis[i]);
				_field_rns[i].invin(_MMi[i]);
				integer tmp= _Mi[i]*(unsigned long)_MMi[i];
				for(size_t j=0;j<_ldm;j++){
					_crt_out[j+i*_ldm]=double(tmp&MASK);
					tmp>>=16;
				}
				double beta=double(1UL<<16);
				double  acc=1;
				for(size_t j=0;j<_ldm;j++){
					_crt_in[j+i*_ldm]=acc;
					_field_rns[i].mulin(acc,beta);
				}
			}
		}

		// Arns must be an array of m*n*_size
		// abs(||A||) <= maxA
		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda,
			  const integer& maxA, bool RNS_MAJOR=false) const
		{
			init(m,n,Arns,rda,A,lda, maxA.bitsize()/16 + (maxA.bitsize()%16?1:0),RNS_MAJOR);
		}

		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const;
		void init_transpose(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const;
		void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false) const;
		void convert_transpose(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false) const;

		
		// reduce entries of Arns to be less than the rns basis elements
		void reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR=false) const;
		

		
	}; // end of struct rns_double

} // end of namespace FFPACK

#include "rns-double.inl"

namespace FFLAS {

	template<>
	inline void fflas_delete (FFPACK::rns_double_elt_ptr A) {FFLAS::fflas_delete( A._ptr);}
	template<>
	inline void fflas_delete (FFPACK::rns_double_elt_cstptr A) {delete[] A._ptr;}

}

#endif // __FFPACK_rns_double_H

