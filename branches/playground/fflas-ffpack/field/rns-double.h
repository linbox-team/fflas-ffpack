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
#include <vector>
using namespace std;

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/field/integer.h"
#include "fflas-ffpack/field/modular-double.h"

// activate only if FFLAS-FFPACK haves multiprecision integer
#ifdef __FFLASFFPACK_HAVE_INTEGER


namespace FFPACK {

	// forward declaration
	struct rns_double_elt_ptr;

	// element of the rns structure (allow virtualization of element from an array of double)
	struct rns_double_elt {
		double    *_ptr;
		size_t  _stride;
		bool     _alloc; // specify wether Element owns its memory; alloc is true only through F.init() and _ptr==NULL (this is to handle Element allocated within a matrix)
		rns_double_elt(): _ptr(NULL), _alloc(false){}
		~rns_double_elt(){ if (_alloc) delete[] _ptr; }

		rns_double_elt(double* p, size_t r, size_t a=false) : _ptr(p), _stride(r), _alloc(a) {}

		inline  rns_double_elt_ptr operator&() ;

		rns_double_elt(const rns_double_elt& x) : _ptr(x._ptr),_stride(x._stride),_alloc(false) {}

	};

	// pointer to element of the rns structure (allow virtualization of element from an array of double)
	struct rns_double_elt_ptr : public rns_double_elt {
		rns_double_elt_ptr(){}
		rns_double_elt_ptr(double* p, size_t r) : rns_double_elt(p,r,false){}

		inline const rns_double_elt& operator*() const {
			return static_cast<const rns_double_elt&>(*this);
		}
		inline  rns_double_elt& operator*()  {
			return static_cast<rns_double_elt&>(*this);
		}
		inline rns_double_elt operator[](size_t i)  const {
			return *((*this)+i);
		}
		inline rns_double_elt_ptr operator+(size_t inc) const{
			return rns_double_elt_ptr(_ptr+inc,_stride);
		}
	};

	rns_double_elt_ptr rns_double_elt::operator&() {
		return 	rns_double_elt_ptr(_ptr,_stride);
	}


	template<>
	rns_double_elt_ptr fflas_const_cast (const rns_double_elt_ptr x){return x;}


	/* Structure that handles rns representation given a bound and bitsize for prime moduli
	 * support sign representation (i.e. the bound must be twice larger then ||A||)
	 */
	struct rns_double {
		typedef Modular<double> ModField;
		vector<double>       _basis; // the rns moduli (mi)
		vector<double>    _invbasis; // the inverse of rns moduli (1/mi)
		vector<ModField> _field_rns; // the associated prime field for each mi
		integer                  _M; // the product of the mi's
		vector<integer>         _Mi; // _M/mi
		vector<double>         _MMi; // (_Mi)^(-1) mod mi
		vector<double>      _crt_in; //  2^(16*j) mod mi
		vector<double>     _crt_out; //  (_Mi._MMi) written in base 2^16
 		size_t                _size; // the size of the rns basis (number of mi's)
		size_t               _pbits; // the size in bit of the mi's
		size_t                 _ldm; // log[2^16](_M)

		typedef double                        BasisElement;
		typedef rns_double_elt                     Element;
		typedef rns_double_elt_ptr             Element_ptr;
		typedef const rns_double_elt_ptr  ConstElement_ptr;

		rns_double(const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
			:  _M(1), _size(0), _pbits(pbits)
		{
			integer::seeding(seed);
			integer prime;
			integer sum=1;
			while (_M < bound*sum) {
				_basis.resize(_size+1);
				do {
					integer::random_exact_2exp(prime, _pbits-1);
					nextprime(prime, prime);
				} while (_M%prime == 0);
				_basis[_size]=prime;
				_size++;
				_M*=prime;
				if (rnsmod) sum+=prime;
			}
			precompute_cst(); 
		}
		
		rns_double(const vector<double>& basis, bool rnsmod=false, long seed=time(NULL))
			:  _basis(basis), _M(1), _size(basis.size()), _pbits(0) 
		{
			for(size_t i=0;i<_size;i++){			
				_M*=_basis[i];
				_pbits=std::max(_pbits, integer(_basis[i]).bitsize());
			}
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

		// Arns must be an array of m*n*_size 
		// abs(||A||) < 2^(16k) 
		void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false) const
		{
			if (k>_ldm)
				FFPACK::Failure(__func__,__FILE__,__LINE__,"rns_struct: init (too large entry)");

			size_t mn=m*n;
			double *A_beta = FFLAS::fflas_new<double >(mn*k);
			const integer* Aiter=A;
			// split A into A_beta according to a Kronecker transform in base 2^16
			for(size_t i=0;i<m;i++)
				for(size_t j=0;j<n;j++){
					size_t idx=j+i*n;
					const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
					const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
					size_t l=0;
					size_t maxs=min(k,(Aiter[j+i*lda].size())<<2);
					if (m0[0]->_mp_size >= 0)
						for (;l<maxs;l++)
							A_beta[l+idx*k]=  m0_ptr[l];
					else
						for (;l<maxs;l++)
							A_beta[l+idx*k]= - double(m0_ptr[l]);
					for (;l<k;l++)
						A_beta[l+idx*k]=  0.;
				}
			if (RNS_MAJOR==false) {
				// Arns = _crt_in x A_beta^T
				cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)_size,(int)mn,(int)k,1.0,_crt_in.data(),(int)_ldm,A_beta,(int)k,0.,Arns,(int)rda);
				// reduce each row i of Arns modulo moduli[i]
				for(size_t i=0;i<_size;i++)
					FFLAS::finit(_field_rns[i],mn,Arns+i*rda,1);
			}
			else {
				// Arns =  A_beta x _crt_in^T
				cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)mn,(int)_size,(int)k,1.0,A_beta,(int)k,_crt_in.data(),(int)_ldm,0.,Arns,(int)_size);
				// reduce each column j of Arns modulo moduli[i]
				for(size_t i=0;i<_size;i++)
					FFLAS::finit(_field_rns[i],mn,Arns+i,_size);
			}
			delete[] A_beta;

#ifdef CHECK_RNS
			bool ok=true;
			for (size_t i=0;i<m;i++)
				for(size_t j=0;j<n;j++)
					for(size_t k=0;k<_size;k++){
						ok= (A[i*lda+j] % (long) _basis[k]) == (long) Arns[i*n+j+k*rda];
						//cout<<A[i*lda+j]<<" mod "<<(long) _basis[k]<<"="<<(long) Arns[i*n+j+k*m*n]<<endl;
					}
			cout<<"RNS finit ... "<<(ok?"OK":"ERROR")<<endl;
#endif
		}


		void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda,
			     const double* Arns, size_t rda, bool RNS_MAJOR=false) const
		{
			integer hM= (_M-1)>>1;
			size_t  mn= m*n;
			double *A_beta= FFLAS::fflas_new<double>(mn*_ldm);

			if (RNS_MAJOR==false)
				// compute A_beta = Ap^T x M_beta
				cblas_dgemm(CblasRowMajor,CblasTrans, CblasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , Arns,(int) rda, _crt_out.data(),(int) _ldm, 0., A_beta,(int)_ldm);
			else // compute A_beta = Ap x M_Beta
				cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)mn, (int)_ldm, (int)_size, 1.0 , Arns, (int)_size, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);

			// compute A using inverse Kronecker transform of A_beta expressed in base 2^log_beta
			integer* Aiter= A;
			size_t k=_ldm;
			size_t k4=((k+3)>>2)+ (((k+3)%4==0)?0:1);
			vector<uint16_t> A0(k4<<2,0),A1(k4<<2,0),A2(k4<<2,0),A3(k4<<2,0);
			integer a0,a1,a2,a3,res;
			mpz_t *m0,*m1,*m2,*m3;
			m0= reinterpret_cast<mpz_t*>(&a0);
			m1= reinterpret_cast<mpz_t*>(&a1);
			m2= reinterpret_cast<mpz_t*>(&a2);
			m3= reinterpret_cast<mpz_t*>(&a3);
			mp_limb_t *m0_d,*m1_d,*m2_d,*m3_d;
			m0_d = m0[0]->_mp_d;
			m1_d = m1[0]->_mp_d;
			m2_d = m2[0]->_mp_d;
			m3_d = m3[0]->_mp_d;
			m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = (int) k4;
			m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size  = m3[0]->_mp_size  = (int) k4;
			for(size_t i=0;i<m;i++)
				for (size_t j=0;j<n;j++){
					size_t idx=i*n+j;
					for (size_t l=0;l<k;l++){
						uint64_t tmp=(uint64_t)A_beta[l+idx*k];
						uint16_t* tptr= reinterpret_cast<uint16_t*>(&tmp);
						A0[l  ]= tptr[0];
						A1[l+1]= tptr[1];
						A2[l+2]= tptr[2];
						A3[l+3]= tptr[3];
					}
					// see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
					m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A0[0]);
					m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A1[0]);
					m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A2[0]);
					m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A3[0]);
					res = a0;res+= a1;res+= a2;res+= a3;
					res%=_M;

					// get the correct result according to the expected sign of A
					if (res>hM)
						res-=_M;
					if (gamma==0)
						Aiter[j+i*lda]=res;
					else
						if (gamma==integer(1))
							Aiter[j+i*lda]+=res;
						else
							if (gamma==integer(-1))
								Aiter[j+i*lda]=res-Aiter[j+i*lda];
							else{
								Aiter[j+i*lda]*=gamma;
								Aiter[j+i*lda]+=res;
							}

				}
			m0[0]->_mp_d = m0_d;
			m1[0]->_mp_d = m1_d;
			m2[0]->_mp_d = m2_d;
			m3[0]->_mp_d = m3_d;
			m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc= m3[0]->_mp_alloc = 1;
			m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size = m3[0]->_mp_size  = 0;

			delete[] A_beta;
#ifdef CHECK_RNS
			bool ok=true;
			for (size_t i=0;i<m;i++)
				for(size_t j=0;j<n;j++)
					for(size_t k=0;k<_size;k++){
						ok= (A[i*lda+j] % (long) _basis[k]) == (long) Arns[i*n+j+k*rda];
						//cout<<A[i*lda+j]<<" mod "<<(long) _basis[k]<<"="<<(long) Arns[i*n+j+k*m*n]<<endl;
					}
			cout<<"RNS finit ... "<<(ok?"OK":"ERROR")<<endl;
#endif

		}
	}; // end of struct rns_double


} // end of namespace FFPACK

namespace FFLAS {
	template<>
	inline void fflas_delete (FFPACK::rns_double_elt_ptr A) {delete[] A._ptr;}

}

#endif
#endif
