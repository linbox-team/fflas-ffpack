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

#include <iterator>     // std::ostream_iterator

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

  /* Structure that handles rns representation with moduli stored in double FP numbers
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
    integer                  _mi_sum; // the product of the mi's

    typedef double                        BasisElement;
    typedef rns_double_elt                     Element;
    typedef rns_double_elt_ptr             Element_ptr;
    typedef rns_double_elt_cstptr     ConstElement_ptr;


    // construct an RNS basis with primes of bit-length (pbits) ensuring to represent integers lying in [0, bound[
    // Rmk: when rnsmod is set to true the RNS basis (m1,m2, ..., mk) satisfies that : m1*m2*...*mk >= bound * (m1+m2+...+mk)
    rns_double(const integer& bound, size_t pbits, bool rnsmod=false, long seed=time(NULL))
      :  _M(1), _size(0), _pbits(pbits), _mi_sum(1)
    {
      integer::seeding(seed);
      Givaro::IntPrimeDom IPD;
      integer prime;
      while (_M < bound*_mi_sum) {
        _basis.resize(_size+1);
        do {
          integer::random_exact_2exp(prime, _pbits-1);
          IPD.nextprimein(prime);
        } while (_M%prime == 0);
        _basis[_size]=prime;
        _size++;
        _M*=prime;
        if (rnsmod) _mi_sum+=prime;
      }
      precompute_cst();
    }

    // construct an RNS basis with (size) primes of bit-length (pbits) 
    rns_double(size_t pbits, size_t size, long seed=time(NULL))
      :  _M(1), _size(size), _pbits(pbits), _mi_sum(1)
    {
      integer::seeding(seed);
      Givaro::IntPrimeDom IPD;
      integer prime;
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

    // construct an RNS basis from a vector of relatively prime numbers (basis)
    template<typename Vect>
    rns_double(const Vect& basis, long seed=time(NULL))
      :  _basis(basis.begin(),basis.end()), _basisMax(basis.size()), _negbasis(basis.size()), _M(1), _size(basis.size()), _pbits(0), _mi_sum(1)
    {
      for(size_t i=0;i<_size;i++){
        _M*=_basis[i];
        _pbits=std::max(_pbits, integer(_basis[i]).bitsize());
      }
      precompute_cst();
    }

    rns_double(const RNSIntegerMod<rns_double>& basis, bool rnsmod=false, long seed=time(NULL)) {

    }

    // can force to reduce integer entries larger than M
    void precompute_cst(size_t K=0){

      // Check that _pbits satisfies log(_M)/16* 2^(_pbits) * 2^16 < 2^53 (Required for correctness)
      // => pbits <= 41 - loglog(_M)
      if ( _pbits > 41 - log(double(_M.bitsize()))/log(2.)){
        std::cout<<"FFLAS Error in rns_double: primes bitsize is too large ... aborting\n";
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
      _crt_in.resize(_size*_ldm);
      _crt_out.resize(_size*_ldm);
      for (size_t i=0;i<_size;i++){
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
        double beta=double(1<<16);
        double  acc=1;
        for(size_t j=0;j<_ldm;j++){
          _crt_in[j+i*_ldm]=acc;
          _field_rns[i].mulin(acc,beta);
        }
      }
    }

    // rda is the distance between two consecutive residues : i.e. A_ij mod mk and A_ij mod m_(k+1) -> it is either 1 (RNS_MAJOR) or mxn (NOT RNS_MAJOR)
    void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false, const FFLAS::FFLAS_TRANSPOSE ta=FFLAS::FflasNoTrans) const;
    void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false, const FFLAS::FFLAS_TRANSPOSE ta=FFLAS::FflasNoTrans) const;

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


  }; // end of struct rns_double



  template<typename RNS>
  class rnsRandIter {
    std::vector<typename RNS::ModField::RandIter> _RNS_rand;
    const RNS& _domain;

  public:
    rnsRandIter(const RNS& R, uint64_t seed=0)
      : _domain(R) {
      for(const auto& F : R._field_rns)
        _RNS_rand.emplace_back(F,seed);
    }

    /** RNS ring Element random assignement.
     * Element is supposed to be initialized
     * @return random ring Element
     */
    typename RNS::Element& random(typename RNS::Element& elt) const {
      auto coefficient(elt._ptr);
      for(auto & iter : _RNS_rand) {
        iter.random( *coefficient );
        coefficient += elt._stride;
      }
      return elt;
    }

    typename RNS::Element& operator()(typename RNS::Element& elt) const {
      return this->random(elt);
    }

    typename RNS::Element operator()() const {
      typename RNS::Element tmp; _domain.init(tmp);
      return this->operator()(tmp);
    }
    typename RNS::Element random() const {
      return this->operator()();
    }

    const RNS& ring() const { return _domain; }

  };

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
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
