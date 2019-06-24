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
        integer                  _mi_sum; // the product of the mi's

        typedef double                        BasisElement;
        typedef rns_double_elt                     Element;
        typedef rns_double_elt_ptr             Element_ptr;
        typedef rns_double_elt_cstptr     ConstElement_ptr;

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
            // std::ostream_iterator<uint64_t> out_it (std::cout,", ");
            // std::cout<<"RNS basis =";
            // std::copy ( _basis.begin(), _basis.end(), out_it );
            // std::cout<<std::endl;
            // std::cout<<"RNS sum Mi ="<<_mi_sum<<"\n";
            precompute_cst();
        }

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

        template<typename Vect>
        rns_double(const Vect& basis, bool rnsmod=false, long seed=time(NULL))
        :  _basis(basis.begin(),basis.end()), _basisMax(basis.size()), _negbasis(basis.size()), _M(1), _size(basis.size()), _pbits(0), _mi_sum(1)
        {
            for(size_t i=0;i<_size;i++){
                //std::cout<<"basis["<<i<<"]="<<_basis[i]<<std::endl;
                _M*=_basis[i];
                _pbits=std::max(_pbits, integer(_basis[i]).bitsize());
            }
            //std::cout<<"M="<<_M<<std::endl;
            precompute_cst();
        }

        rns_double(const RNSIntegerMod<rns_double>& basis, bool rnsmod=false, long seed=time(NULL)) {

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
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                for(;l<maxs;l++)
                    _crt_out[l+i*_ldm]=m0_ptr[l];
#else
                for(;l<maxs;l++)
                    _crt_out[l+i*_ldm]=m0_ptr[l ^ ((sizeof(mp_limb_t)/2U) - 1U)];
#endif
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
        std::vector<double>      _crt_in; //  2^(16*j) mod mi
        std::vector<double>     _crt_out; //  (_Mi._MMi) written in base 2^16
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
            _ldm = (_M.bitsize()/16) + ((_M.bitsize()%16)?1:0) ;
            _invbasis.resize(_size);
            _basisMax.resize(_size);
            _negbasis.resize(_size);
            _field_rns.resize(_size);
            _Mi.resize(_size);
            _MMi.resize(_size);
            _crt_in.resize(_size*_ldm);
            _crt_out.resize(_size*_ldm);
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
                for(size_t j=0;j<_ldm;j++){
                    _crt_out[j+i*_ldm]=double(tmp&MASK);
                    tmp>>=16;
                }
                double beta=double(1<<16);
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
            init(m*n,Arns,A,lda);
        }

        void init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR=false){
            init(m*n,Arns,A,lda);
        }
        void convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR=false){
            convert(m*n, A, Arns);
        }
        void init(size_t m, double* Arns, const integer* A, size_t lda) const;
        void convert(size_t m, integer *A, const double *Arns) const;

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)

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

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

        // reduce entries of Arns to be less than the rns basis elements
        void reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR=false) const;



    }; // end of struct rns_double_extended


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
