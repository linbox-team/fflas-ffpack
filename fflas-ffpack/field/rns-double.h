/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
#include <givaro/givintprime.h>
#include "givaro/modular-extended.h"
#include <recint/ruint.h>
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
		std::vector<double, AlignedAllocator<double, Alignment::CACHE_LINE>>       _invbasis; // the inverse of rns moduli (1/mi)
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

		rns_double(size_t pbits, size_t size, long seed=time(NULL))
		:  _M(1), _size(size), _pbits(pbits)
		{
			integer::seeding(seed);
                        Givaro::IntPrimeDom IPD;
			integer prime;
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
		rns_double(const Vect& basis, bool rnsmod=false, long seed=time(NULL))
			:  _basis(basis.begin(),basis.end()), _basisMax(basis.size()), _negbasis(basis.size()), _M(1), _size(basis.size()), _pbits(0)
		{
			for(size_t i=0;i<_size;i++){
				//std::cout<<"basis["<<i<<"]="<<_basis[i]<<std::endl;
				_M*=_basis[i];
				_pbits=std::max(_pbits, integer(_basis[i]).bitsize());
			}
			//std::cout<<"M="<<_M<<std::endl;
			precompute_cst();
		}

		// can force to reduce integer entries larger than M
		void precompute_cst(size_t K=0){
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
			_crt_in.resize(_size*_ldm);
			_crt_out.resize(_size*_ldm);
			//const unsigned int MASK=0xFFFF;
			//Givaro::Timer chrono;
			//double t1=0.,t2=0.,t3=0.;
			
			for (size_t i=0;i<_size;i++){
				//chrono.start();
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
				//chrono.stop();
				//t1+=chrono.usertime();
				//chrono.start();
				/*
				for(size_t j=0;j<_ldm;j++){
					_crt_out[j+i*_ldm]=double(tmp[0]&MASK);
					tmp>>=16; // Bad idea -> too slow (must get the lowest limb of the integer)
					
				}
				*/
				size_t l=0;
				for(;l<maxs;l++)
					_crt_out[l+i*_ldm]=m0_ptr[l];
				for(;l<_ldm;l++)
					_crt_out[l+i*_ldm]=0.;;
				// chrono.stop();
				// t2+=chrono.usertime();
				// chrono.start();			       				
				double beta=double(1<<16);
				double  acc=1;	       
				for(size_t j=0;j<_ldm;j++){
					_crt_in[j+i*_ldm]=acc;
					_field_rns[i].mulin(acc,beta);					
				}
				// chrono.stop();
				// t3+=chrono.usertime();
			
			}
			// std::cout<<"t1="<<t1<<std::endl;
			// std::cout<<"t2="<<t2<<std::endl;
			// std::cout<<"t3="<<t3<<std::endl;
		 }

		// Arns must be an array of m*n*_size
		// abs(||A||) <= maxA
		template<typename T>
		void init(size_t m, size_t n, double* Arns, size_t rda, const T* A, size_t lda,
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

		template<size_t K>
		void init(size_t m, size_t n, double* Arns, size_t rda, const RecInt::ruint<K>* A, size_t lda, size_t k, bool RNS_MAJOR=false) const;
		template<size_t K>
		void convert(size_t m, size_t n, integer gamma, RecInt::ruint<K>* A, size_t lda, const double* Arns, size_t rda, integer p=0,bool RNS_MAJOR=false) const;

		
	}; // end of struct rns_double
	
	/* Structure that handles rns representation given a bound and bitsize for prime moduli, allow large moduli
	 * support sign representation (i.e. the bound must be twice larger then ||A||)
	 */
	struct rns_double_extended {
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

		rns_double_extended(const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
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

		rns_double_extended(size_t pbits, size_t size, long seed=time(NULL))
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
		rns_double_extended(const Vect& basis, bool rnsmod=false, long seed=time(NULL))
			:  _basis(basis.begin(),basis.end()), _basisMax(basis.size()), _negbasis(basis.size()), _M(1), _size(basis.size()), _pbits(0)
		{
			for(size_t i=0;i<_size;i++){
				//std::cout<<"basis["<<i<<"]="<<_basis[i]<<std::endl;
				_M*=_basis[i];
				_pbits=std::max(_pbits, integer(_basis[i]).bitsize());
			}
			//std::cout<<"M="<<_M<<std::endl;
			precompute_cst();
		}


		void precompute_cst(){
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
					_crt_out[3][j+i*_ldm]=a0-a2;
					_crt_out[4][j+i*_ldm]=a0+a1+a2;
					_crt_out[5][j+i*_ldm]=a0+a1-a2;					
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
			std::cout<<"basis:= [";
			for(size_t i=0;i<_size;i++)
				std::cout<<(int64_t)_basis[i]<<(i!=_size-1?",":"];\n");
		
			for(size_t l=0;l<6;l++){
				std::cout<<"crtin"<<l<<":=Matrix("<<_size<<","<<_ldm<<",[";
				for(size_t i=0;i<_size;i++){
					std::cout<<"[";
					for(size_t j=0;j<_ldm;j++)
						std::cout<<(int64_t)_crt_in[l][j+i*_ldm]<<(j!=_ldm-1?",":(i!=_ldm-1?"],":"]"));
				}
				std::cout<<"]);\n";
			}
#endif
			

		}

		// Arns must be an array of m*n*_size
		// abs(||A||) <= maxA
		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda,
			  const integer& maxA, bool RNS_MAJOR=false) const
		{
			init(m,n,Arns,rda,A,lda, uint64_t(maxA.bitsize()/48 + (maxA.bitsize()%48?1:0)),RNS_MAJOR);
		}
		

		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const
		{
		    //init(m*n,Arns,A,lda);
			if (k>_ldm){
				FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_struct: init (too large entry)");
				std::cerr<<"k="<<k<<" _ldm="<<_ldm<<std::endl;
		}
		size_t mn=m*n;
		size_t mnk=mn*k;
		double *A_beta = FFLAS::fflas_new<double >(mnk*6);
		double *A_beta0=A_beta;
		double *A_beta1=A_beta+1*mnk;
		double *A_beta2=A_beta+2*mnk;
		double *A_beta3=A_beta+3*mnk;
		double *A_beta4=A_beta+4*mnk;
		double *A_beta5=A_beta+5*mnk;
		size_t mnsize=mn*_size;
		double *A_rns_tmp = FFLAS::fflas_new<double >(mnsize*6);
		double *A_rns0=A_rns_tmp;
		double *A_rns1=A_rns_tmp+1*mnsize;
		double *A_rns2=A_rns_tmp+2*mnsize;
		double *A_rns3=A_rns_tmp+3*mnsize;
		double *A_rns4=A_rns_tmp+4*mnsize;
		double *A_rns5=A_rns_tmp+5*mnsize;

		
		const integer* Aiter=A;
		// split A into A_beta according to a Kronecker transform in base 2^48
		Givaro::Timer tkr; tkr.start();
		PARFOR1D(i,m,SPLITTER(NUM_THREADS),
				 for(size_t j=0;j<n;j++){
					 size_t idx=j+i*n;
					  const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
					  const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
					  size_t l=0,h=0;
					  // size in base 2^16
					  size_t k16=((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2);
					  size_t maxs=std::min(k,((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3);// to ensure 32 bits portability
					  double a0,a1,a2;
					  if (m0[0]->_mp_size >= 0)
						  for (;l<maxs;l++){
							  a0= m0_ptr[h++];//std::cout<<"a0="<<(uint64_t)a0<<std::endl;
							  a1= m0_ptr[h++];//std::cout<<"a1="<<(uint64_t)a1<<std::endl;
							  a2= m0_ptr[h++];//std::cout<<"a2="<<(uint64_t)a2<<std::endl;
							  A_beta0[l+idx*k]= a0;
							  A_beta1[l+idx*k]= a1;
							  A_beta2[l+idx*k]= a2;
							  A_beta3[l+idx*k]= a0+a1;
							  A_beta4[l+idx*k]= a1+a2;
							  A_beta5[l+idx*k]= a0+a1+a2;							  							  
						  }
					  else
						  for (;l<maxs;l++){
							  a0= -double(m0_ptr[h++]);
							  a1= -double(m0_ptr[h++]);
							  a2= -double(m0_ptr[h++]);
							  A_beta0[l+idx*k]= a0;
							  A_beta1[l+idx*k]= a1;
							  A_beta2[l+idx*k]= a2;
							  A_beta3[l+idx*k]= a0+a1;
							  A_beta4[l+idx*k]= a1+a2;
							  A_beta5[l+idx*k]= a0+a1+a2;							  							  
						  }
					  for (;l<k;l++){
						  a0= (h<k16)?m0_ptr[h++]:0.;
						  a1= (h<k16)?m0_ptr[h++]:0.;
						  a2= (h<k16)?m0_ptr[h++]:0.;
						  A_beta0[l+idx*k]= a0;
						  A_beta1[l+idx*k]= a1;
						  A_beta2[l+idx*k]= a2;
						  A_beta3[l+idx*k]= a0+a1;
						  A_beta4[l+idx*k]= a1+a2;
						  A_beta5[l+idx*k]= a0+a1+a2;							  							  							  
					  }					
				 }
				 );

#ifdef CHECK_RNS
		for (size_t i=0;i<m;i++)
			for (size_t j=0;j<n;j++){
				//std::cout<<"A="<<A[i*lda+j]<<std::endl;
				
				
				integer tmp=0;
				int idx=j+i*n;
				for(int l=k-1;l>=0;l--){
					int64_t limb,c0,c1,c2,c3,c4;
					c0= A_beta0[l+idx*k]; // 1
					c1= A_beta3[l+idx*k] - A_beta1[l+idx*k] - A_beta0[l+idx*k] ; // X
					c2= A_beta5[l+idx*k] - A_beta3[l+idx*k] - A_beta4[l+idx*k] + 2*A_beta1[l+idx*k]; // X^2
					c3= A_beta4[l+idx*k] - A_beta1[l+idx*k] - A_beta2[l+idx*k]; // X^3
					c4= A_beta2[l+idx*k]; // X^4
					if (c1!=0. || c3!=0.) std::cout<<"RNS EXTENDED: error in splitting entries (linear form)\n";
					limb= c0+(c2<<16)+(c4<<32);
					tmp=(tmp<<48)+limb;
				}
				if (tmp!=A[i*lda+j]){
					std::cout<<"RNS EXTENDED: error in splitting entries (16-adic) --> "<<tmp<<" != "<<A[i*lda+j]<<std::endl;;
					std::cout<<"MAX="<<std::min(k,((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3)<<std::endl;
					std::cout<<"MAX1="<<k<<std::endl;
					std::cout<<"MAX2="<<((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3<<std::endl;
					std::cout<<"bisize: "<<Aiter[j+i*lda].bitsize()<<std::endl;

					tmp=A[i*lda+j];
					for(int l=0;l<k;l++){
						std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFF))             << " != "<<int64_t(A_beta0[l+idx*k])<<std::endl;
						std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFFFFFF)>>16)     << " != "<<int64_t(A_beta1[l+idx*k])<<std::endl;
						std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFFFFFFFFFF)>>32) << " != "<<int64_t(A_beta2[l+idx*k])<<std::endl;
						tmp>>=48;
					}
				}
			}													
#endif
			tkr.stop();
			std::cerr<<"RNS double ext - Kronecker : "<<tkr<<std::endl;
			if (RNS_MAJOR==false) {
				// A_rns = _crt_in x A_beta^T
				Givaro::Timer tfgemm; tfgemm.start();
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[0].data(),_ldm,A_beta0,k,0.,A_rns0, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[1].data(),_ldm,A_beta1,k,0.,A_rns1, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[2].data(),_ldm,A_beta2,k,0.,A_rns2, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[3].data(),_ldm,A_beta3,k,0.,A_rns3, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[4].data(),_ldm,A_beta4,k,0.,A_rns4, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
				FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[5].data(),_ldm,A_beta5,k,0.,A_rns5, mn,
							  FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());			
				tfgemm.stop();
				std::cerr<<"RNS double ext - fgemm : "<<tfgemm<<std::endl;
			}
			else {
				// Arns =  A_beta x _crt_in^T
				std::cerr<<"NOT YET IMPLEMENTED .... aborting\n"; std::terminate();
			}
			Givaro::Timer tred; tred.start();

#ifdef RNS_DEBUG			
			for(size_t l=0;l<6;l++){
				std::cout<<"Abeta"<<l<<":=Matrix("<<_size<<","<<mn<<",[";
				for(size_t i=0;i<k;i++){
					std::cout<<"[";
					for(size_t j=0;j<mn;j++)
						std::cout<<(int64_t)A_beta[l*mnk+j+i*mn]<<(j!=mn-1?",":(i!=k-1?"],":"]"));
				}
				std::cout<<"]);\n";
			}

			for(size_t l=0;l<6;l++){
				std::cout<<"Arns"<<l<<":=Matrix("<<_size<<","<<mn<<",[";
				for(size_t i=0;i<_size;i++){
					std::cout<<"[";
					for(size_t j=0;j<mn;j++)
						std::cout<<(int64_t)A_rns_tmp[l*mnsize+j+i*mn]<<(j!=mn-1?",":(i!=_size-1?"],":"]"));
				}
				std::cout<<"]);\n";
			}
#endif
			
			
			double c,c0,c1,c2,c3,c4;
			PARFOR1D(i,_size,SPLITTER(NUM_THREADS),
					 double twosixty_mod_mi;
					 _field_rns[i].init(twosixty_mod_mi,(1<<16));
					 for(size_t j=0;j<mn;j++)
						 {
							 size_t h=j+i*mn;
							 c0= A_rns0[h]; // 1
							 c1= A_rns3[h] - A_rns1[h] - A_rns0[h] ; // X
							 c2= A_rns5[h] - A_rns3[h] - A_rns4[h] + 2*A_rns1[h]; // X^2
							 c3= A_rns4[h] - A_rns1[h] - A_rns2[h]; // X^3
							 c4= A_rns2[h]; // X^4
#ifdef RNS_DEBUG
							 std::cout<<(int64_t)c0
									  <<"+2^16*"<<(int64_t)c1
									  <<"+2^32*"<<(int64_t)c2
									  <<"+2^48*"<<(int64_t)c3
									  <<"+2^64*"<<(int64_t)c4<<";\n";
#endif				 
							 // compute c=c0+c1.2^16+c2.2^32+c3.2^48+c4.2^64 mod mi
							 _field_rns[i].axpy(c,c4,twosixty_mod_mi,c3);
							 _field_rns[i].axpy(c,c ,twosixty_mod_mi,c2);
							 _field_rns[i].axpy(c,c ,twosixty_mod_mi,c1);
							 _field_rns[i].axpy(c,c ,twosixty_mod_mi,c0);
							 Arns[j+i*rda]= c+(c>=0?0:_basis[i]);
						 }
					 );
			//reduce(mn,Arns,rda,RNS_MAJOR);
			
			tred.stop();

			std::cerr<<"RNS double ext - Reduce : "<<tred<<std::endl;
	
		FFLAS::fflas_delete(A_beta);
		FFLAS::fflas_delete(A_rns_tmp);

#ifdef CHECK_RNS
		bool ok=true;
		for (size_t i=0;i<m;i++)
			for(size_t j=0;j<n;j++)
				for(size_t k=0;k<_size;k++){
					ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0)) == (int64_t) Arns[i*n+j+k*rda]);
					if (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
					    != (int64_t) Arns[i*n+j+k*rda])
					{
						std::cout<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
								 <<" != "
								 <<(int64_t) Arns[i*n+j+k*rda]
								 <<" --> "<<A[i*lda+j]<<"  mod ("<<(uint64_t)_basis[k]<<")"<<std::endl;
					}
					else{
						// std::cout<<"OK -> "<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
						// 		 <<" == "
						// 		 <<(int64_t) Arns[i*n+j+k*rda]
						// 		 <<" --> "<<A[i*lda+j]<<"  mod ("<<(uint64_t)_basis[k]<<")"<<std::endl;
						
					}
				}
		std::cout<<"RNS EXTENDED freduce ... "<<(ok?"OK":"ERROR")<<std::endl;
		if (!ok) std::terminate();
#endif



			
		}






		
		void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false) const {
		  convert(m*n, A, Arns);
		}
		void init(size_t m, double* Arns, const integer* A, size_t lda) const;
		void convert(size_t m, integer *A, const double *Arns) const;
		
#if defined(__FFLASFFPACK_USE_SIMD)
		
		template<class SimdT>
		inline void splitSimd(const SimdT x, SimdT & x_h, SimdT & x_l) const {
			using simd = Simd<double>;
			using vect_t = typename simd::vect_t;
			vect_t vc = simd::set1((double)((1 << 27)+1));
			vect_t tmp = simd::mul(vc, x);
			x_h = simd::add(tmp, simd::sub(x, tmp));
			x_l = simd::sub(x, x_h);
		}

		template<class SimdT>
		inline void multSimd(const SimdT va, const SimdT vb, SimdT & vs, SimdT & vt) const{
			using simd = Simd<double>;
			using vect_t = typename simd::vect_t;
			vect_t vah, val, vbh, vbl;
			vs = simd::mul(va, vb);
//#ifdef __FMA__
			vt = simd::fnmadd(va, vb, vs);
//#else
			splitSimd(va, vah, val);
			splitSimd(vb, vbh, vbl);		
			vt = simd::add(simd::add(simd::sub(simd::mul(vah, vbh), vs), simd::mul(vah, vbl)), simd::add(simd::mul(val, vbh), simd::mul(val, vbl)));
//#endif
		}
		
		template<class SimdT>
		inline SimdT modSimd(const SimdT a, const SimdT p, const SimdT ip, const SimdT np) const{
		  using simd = Simd<double>;
		  using vect_t = typename simd::vect_t;
		  vect_t pqh, pql, abl, abh;
		  vect_t q = simd::floor(simd::mul(a, ip));
		  multSimd(p, q, pqh, pql);
		  vect_t r = simd::add(simd::sub(a, pqh), pql);
		  abh = simd::greater_eq(r, p);
		  abl = simd::lesser(r, simd::zero());
		  abh = simd::vand(abh, np);
		  abl = simd::vand(abl, p);
		  abh = simd::vor(abh, abl);
		  return r = simd::add(r, abh);
		}
		
#endif // __FFLASFFPACK_USE_SIMD
		
		// reduce entries of Arns to be less than the rns basis elements
		void reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR=false) const;
		

		
	}; // end of struct rns_double_extended

} // end of namespace FFPACK

#include "rns-double.inl"
#include "rns-double-recint.inl"
namespace FFLAS {

	template<>
	inline void fflas_delete (FFPACK::rns_double_elt_ptr A) {FFLAS::fflas_delete( A._ptr);}
	template<>
	inline void fflas_delete (FFPACK::rns_double_elt_cstptr A) {delete[] A._ptr;}

}

#endif // __FFPACK_rns_double_H

