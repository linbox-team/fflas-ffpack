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

/*! @file field/rns-double-extended.h
 * @ingroup field
 * @brief  rns structure with double support (full precision)
 */

#ifndef __FFPACK_rns_double_extended_H
#define __FFPACK_rns_double_extended_H

// Bigger multiple of s lesser or equal than x, s must be a power of two
#ifndef ROUND_DOWN
#define ROUND_DOWN(x, s) ((x) & ~((s)-1))
#endif

#include <vector>
#include <givaro/modular-floating.h>
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include "givaro/modular-extended.h"
#include <recint/ruint.h>
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include "fflas-ffpack/field/rns-double-elt.h"

namespace FFPACK {

  // forward declaration of different variants
  struct  rns_double_extended_V1 ;
  struct  rns_double_extended_V2 ;

  // set default value to one variant
  using rns_double_extended = rns_double_extended_V2;


#define RALIGN setw(40)<<right
	/* Structure that handles rns representation given a bound and bitsize for prime moduli, allow large moduli
	 * support sign representation (i.e. the bound must be twice larger then ||A||)
	 */
	struct rns_double_extended_V1 {
		typedef Givaro::Integer integer;
		typedef Givaro::ModularExtended<double> ModField;
		
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basis; // the rns moduli (mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basisMax; // (mi-1)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _negbasis; // (-mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _invbasis; // the inverse of rns moduli (1/mi)
		std::vector<ModField> _field_rns; // the associated prime field for each mi
		integer                  _M; // the product of the mi's
		std::vector<integer>         _Mi; // _M/mi
		std::vector<double>         _MMi; // (_Mi)^(-1) mod mi
		std::vector<double>      _crt_in[6]; //  2^(16*j) mod mi 
		std::vector<double>     _crt_out[6]; //  (_Mi._MMi) written in base 2^16
		size_t                _size; // the size of the rns basis (number of mi's)
		size_t               _pbits; // the size in bit of the mi's
		size_t                 _ldm; // log[2^16](_M)

		typedef double                        BasisElement;
		typedef rns_double_elt                     Element;
		typedef rns_double_elt_ptr             Element_ptr;
		typedef rns_double_elt_cstptr     ConstElement_ptr;

		rns_double_extended_V1(const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
			:  _M(1), _size(0), _pbits(pbits)
		{
			integer::seeding(seed);
			integer prime; Givaro::IntPrimeDom IPD;
			integer sum=1;
			while (_M < bound*sum) {
				_basis.resize(_size+1);
				do {
					integer::random_exact_2exp(prime, _pbits-1);
					IPD.nextprimein(prime);
				} while (_M%prime == 0);
				_basis[_size]=prime;
				_size++;
				_M*=prime;
				if (rnsmod) sum+=prime;
			}
			precompute_cst();
		}

		rns_double_extended_V1(size_t pbits, size_t size, long seed=time(NULL))
			:  _M(1), _size(size), _pbits(pbits)
		{
			integer::seeding(seed);
			integer prime; Givaro::IntPrimeDom IPD;
			integer sum=1;
			_basis.resize(size);
			_negbasis.resize(size);
			_basisMax.resize(size);
			for(size_t i = 0 ; i < _size ; ++i){
				integer::random_exact_2exp(prime, _pbits-1);
				IPD.nextprimein(prime);
				_basis[i]=prime;
				_basisMax[i] = prime-1;
				_negbasis[i] = 0-prime;
				_M*=prime;
			}
			precompute_cst();
		}

		template<typename Vect>
		rns_double_extended_V1(const Vect& basis, bool rnsmod=false, long seed=time(NULL))
			:  _basis(basis.begin(),basis.end()), _basisMax(basis.size()), _negbasis(basis.size()), _M(1), _size(basis.size()), _pbits(0)
		{
			for(size_t i=0;i<_size;i++){
				_M*=_basis[i];
				_pbits=std::max(_pbits, integer(_basis[i]).bitsize());
			}
			precompute_cst();
		}


		void precompute_cst(){

      // Check that _pbits <= 48 
      if ( _pbits > 48 ){
        std::cout<<"FFLAS Error in rns_double_extended: primes bitsize "<<_pbits<<" is too large ... aborting\n";
        std::terminate();
      }

			_ldm = (_M.bitsize()/48) + ((_M.bitsize()%48)?1:0) ;
			_invbasis.resize(_size);
			_basisMax.resize(_size);
			_negbasis.resize(_size);
			_field_rns.resize(_size);
			_Mi.resize(_size);
			_MMi.resize(_size);
			for (size_t i=0;i<6;i++){
				_crt_in[i] .resize(_size*_ldm);
				_crt_out[i].resize(_size*_ldm);
			}
			const unsigned int MASK=0xFFFF;
			for (size_t i=0;i<_size;i++){
				_invbasis[i]  = 1./_basis[i];
				_basisMax[i] = _basis[i]-1;
				_negbasis[i] = 0-_basis[i];
				_field_rns[i] = ModField(_basis[i]);
				_Mi[i]        = _M/(uint64_t)_basis[i];
				_field_rns[i].init(_MMi[i], _Mi[i] % (double)_basis[i]);
				_field_rns[i].invin(_MMi[i]);
				integer tmp= _Mi[i]*(uint64_t)_MMi[i];
				double a0,a1,a2;
				for(size_t j=0;j<_ldm;j++){
					uint64_t limb= tmp[0];
					a0= double(limb&MASK);
					a1= double((limb>>16)&MASK);
					a2= double((limb>>32)&MASK);					
					_crt_out[0][j+i*_ldm]=a0;
					_crt_out[1][j+i*_ldm]=a1;
					_crt_out[2][j+i*_ldm]=a2;
					_crt_out[3][j+i*_ldm]=a0+a1;
					_crt_out[4][j+i*_ldm]=a1+a2;
					_crt_out[5][j+i*_ldm]=a0+a1+a2;					
					tmp>>=48;
				}
				double beta=double(1UL<<48);
				double  acc=1;
				for(size_t j=0;j<_ldm;j++){
					uint64_t limb= acc;
					a0= double(limb&MASK);
					a1= double((limb>>16)&MASK);
					a2= double((limb>>32)&MASK);					
					_crt_in[0][j+i*_ldm]=a0;
					_crt_in[1][j+i*_ldm]=a1;
					_crt_in[2][j+i*_ldm]=a2;
					_crt_in[3][j+i*_ldm]=a0+a1;
					_crt_in[4][j+i*_ldm]=a1+a2;
					_crt_in[5][j+i*_ldm]=a0+a1+a2;					
					_field_rns[i].mulin(acc,beta);
				}
			}
#ifdef RNS_DEBUG
			std::cout<<"basis= [";
			for(size_t i=0;i<_size;i++)
				std::cout<<(int64_t)_basis[i]<<(i!=_size-1?",":"]\n");
      Givaro::ModularExtended<double> ZZ(2UL<<48);;
			for(size_t l=0;l<6;l++){
        std::cout<<"crtin"<<l<<"=\n";
        FFLAS::WriteMatrix(std::cout, ZZ, _size, _ldm, _crt_in[l].data(), _ldm);
        std::cout<<std::endl;				
			}
			for(size_t l=0;l<6;l++){
				std::cout<<"crtout"<<l<<"=\n";
        FFLAS::WriteMatrix(std::cout, ZZ, _size, _ldm, _crt_out[l].data(), _ldm);
        std::cout<<std::endl;				
			}
#endif
		}

		// Arns must be an array of m*n*_size
		// abs(||A||) <= maxA
		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, const integer& maxA, bool RNS_MAJOR=false) const
		{
			init(m,n,Arns,rda,A,lda, uint64_t(maxA.bitsize()/48 + (maxA.bitsize()%48?1:0)),RNS_MAJOR);
		}
	  
    void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const;
		void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false) const;

		// reduce entries of Arns to be less than the rns basis elements
		void reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR=false) const;

	}; // end of struct rns_double_extended





	struct rns_double_extended_V2 {
		typedef Givaro::Integer integer;
		typedef Givaro::ModularExtended<double> ModField;
		
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basis; // the rns moduli (mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _basisMax; // (mi-1)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _negbasis; // (-mi)
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _invbasis; // the inverse of rns moduli (1/mi)
		std::vector<ModField> _field_rns; // the associated prime field for each mi
		integer                  _M; // the product of the mi's
		std::vector<integer>         _Mi; // _M/mi
		std::vector<double>         _MMi; // (_Mi)^(-1) mod mi
		std::vector<double>      _crt_in; //  2^(16*j) mod mi -> each entry are splitted into two line (a0 + 2^27 a1)
		std::vector<double>     _crt_out; //  (_Mi._MMi) written in base 2^16
		size_t                _size; // the size of the rns basis (number of mi's)
		size_t               _pbits; // the size in bit of the mi's
		size_t                 _ldm; // log[2^16](_M)
		uint64_t    _shift;

		typedef double                        BasisElement;
		typedef rns_double_elt                     Element;
		typedef rns_double_elt_ptr             Element_ptr;
		typedef rns_double_elt_cstptr     ConstElement_ptr;

		rns_double_extended_V2 (const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
			:  _M(1), _size(0), _pbits(pbits)
		{
			if (pbits > 52){
				std::cerr<<"FFLAS- RNS EXTENDED ERROR: prime bitsize above 52 bits ... log(p)="<<pbits<<std::endl;
				std::terminate();
			}
			_shift=pbits/2;

			if( (bound.bitsize()/uint64_t(pbits)) * (integer(1)<<_shift) >= (integer(1)<<37) ){
				std::cerr<<"FFLAS- RNS EXTENDED ERROR: prime bitsize ("<<pbits<<") too large to handle RNS basis with "<<bound.bitsize()<<" bits"<<std::endl;
				std::terminate();				
      }


			integer::seeding(seed);
			Givaro::IntPrimeDom IPD;
			integer prime;
			integer sum=1;
			while (_M < bound*sum) {
				_basis.resize(_size+1);
				do {
					integer::random_exact_2exp(prime, _pbits-1);
					IPD.nextprimein(prime);
				} while (_M%prime == 0);
				_basis[_size]=prime;
				_size++;
				_M*=prime;
				if (rnsmod) sum+=prime;
			}
			precompute_cst();
		}
		// can force to reduce integer entries larger than M
		void precompute_cst(size_t K=0){

      // Check that _pbits satisfies log(_M)/16* 2^(_pbits/2) * 2^16 < 2^53 (Required for correctness)
      // => pbits <= 2*(41 - loglog(_M))
      if ( _pbits > 2*(41 - log(double(_M.bitsize()))/log(2.))){
        std::cout<<"FFLAS Error in rns_double_extended: primes bitsize "<<_pbits<<" is too large ... aborting\n";
        std::terminate();
      }

			if (K!=0)
				_ldm=K;
			else
				_ldm = (_M.bitsize()/16) + ((_M.bitsize()%16)?1:0) ;
			_invbasis.resize(_size);
			_field_rns.resize(_size);
			_Mi.resize(_size);
			_MMi.resize(_size);
			_basisMax.resize(_size);
			_negbasis.resize(_size);			
			_crt_in.resize(2*_size*_ldm);
			_crt_out.resize(_size*_ldm);
			//const unsigned int MASK=0xFFFF;
#ifdef BENCH_RNS_PRECOMP
			Givaro::Timer chrono;
			double t1=0.,t2=0.,t3=0.;
#endif
			for (size_t i=0;i<_size;i++){
#ifdef BENCH_RNS_PRECOMP
				chrono.start();
#endif
				_invbasis[i]  = 1./_basis[i];
				_basisMax[i] = _basis[i]-1;
				_negbasis[i] = 0-_basis[i];
				_field_rns[i] = ModField(_basis[i]);
				_Mi[i]        = _M/(uint64_t)_basis[i];
				_field_rns[i].init(_MMi[i], _Mi[i] % (double)_basis[i]);				
				_field_rns[i].invin(_MMi[i]);
				integer tmp= _Mi[i]*(uint64_t)_MMi[i];
				const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(&tmp);
				const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
				size_t maxs=std::min(_ldm,(tmp.size())*sizeof(mp_limb_t)/2);// to ensure 32 bits portability
#ifdef BENCH_RNS_PRECOMP
				chrono.stop();
				t1+=chrono.usertime();
				chrono.start();
#endif 
				size_t l=0;
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
        for(;l<maxs;l++)
          _crt_out[l+i*_ldm]=m0_ptr[l];
#else
        for(;l<maxs;l++)
          _crt_out[l+i*_ldm]=m0_ptr[l ^ ((sizeof(mp_limb_t)/2U) - 1U)];
#endif
        for(;l<_ldm;l++)
					_crt_out[l+i*_ldm]=0.;;
#ifdef BENCH_RNS_PRECOMP
				chrono.stop();
        t2+=chrono.usertime();
        chrono.start();
#endif
        double beta=double(1<<16);
        double  acc=1;

        for(size_t j=0;j<_ldm;j++){
          uint64_t acci= (uint64_t)acc;
          _crt_in[j+i*_ldm]=acci  & ((1<<_shift)-1);
          _crt_in[j+(i+_size)*_ldm]=(acci >> _shift);
          _field_rns[i].mulin(acc,beta);					

        }
#ifdef BENCH_RNS_PRECOMP
				chrono.stop();
				t3+=chrono.usertime();
#endif
			}
#ifdef BENCH_RNS_PRECOMP
			std::cout<<"RNS precomp t1="<<t1<<std::endl;
			std::cout<<"RNS precomp t2="<<t2<<std::endl;
			std::cout<<"RNS precomp t3="<<t3<<std::endl;
#endif


#ifdef RNS_DEBUG
			std::cout<<"RNS double ext - basis:= [";
			for(size_t i=0;i<_size;i++)
				std::cout<<(int64_t)_basis[i]<<(i!=_size-1?",":"];\n");
			Givaro::ModularExtended<double> ZZ(2UL<<48);;
			std::cout<<"CRTmat:=";
      //write_field(ZZ,std::cout,_crt_out.data(), _size, _ldm,_ldm,true);
      FFLAS::WriteMatrix(std::cout, ZZ, _size, _ldm, _crt_out.data(), _ldm);
#endif			
		}

    // Arns must be an array of m*n*_size
    // abs(||A||) <= maxA    



    void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false, const FFLAS::FFLAS_TRANSPOSE ta=FFLAS::FflasNoTrans) const;
    void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false, const FFLAS::FFLAS_TRANSPOSE ta=FFLAS::FflasNoTrans) const;
    void convert_bis(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false, const FFLAS::FFLAS_TRANSPOSE ta=FFLAS::FflasNoTrans) const;

    // Arns must be an array of m*n*_size with  abs(||A||) <= maxA    
    void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, const integer& maxA, bool RNS_MAJOR=false) const{
      init(m,n,Arns,rda,A,lda, maxA.bitsize()/16 + (maxA.bitsize()%16?1:0),RNS_MAJOR, FFLAS::FflasNoTrans);
    }
    void init_transpose(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, const integer& maxA, bool RNS_MAJOR=false) const{
      init(m,n,Arns,rda,A,lda, maxA.bitsize()/16 + (maxA.bitsize()%16?1:0),RNS_MAJOR, FFLAS::FflasTrans);
    }
    void init_transpose(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const{
      init(m,n,Arns,rda,A,lda,k,RNS_MAJOR, FFLAS::FflasTrans);
    }
    void convert_transpose(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false) const {
      convert(m,n,gamma,A,lda,Arns,rda,RNS_MAJOR,FFLAS::FflasTrans);
    }



    // reduce entries of Arns to be less than the rns basis elements
    void reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR=false) const;

    template<size_t K>
    void init(size_t m, size_t n, double* Arns, size_t rda, const RecInt::ruint<K>* A, size_t lda, size_t k, bool RNS_MAJOR=false) const;
    template<size_t K>
    void convert(size_t m, size_t n, integer gamma, RecInt::ruint<K>* A, size_t lda, const double* Arns, size_t rda, integer p=0,bool RNS_MAJOR=false) const;
    
	};


   
	
} // end of namespace FFPACK

#include "rns-double-extended.inl"

#endif //__FFPACK_rns_double_extended_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
