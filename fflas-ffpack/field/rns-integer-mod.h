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
/*! @file field/rns-integer-mod.h
 * @ingroup field
 * @brief  representation of <code>Z/pZ</code> using RNS representation (note: fixed precision)
 */


#ifndef __FFPACK_rns_integer_mod_H
#define __FFPACK_rns_integer_mod_H

#include <vector>
#include <cmath>

#include <recint/recint.h>
#include <givaro/modular-integer.h>
#include <givaro/givinteger.h>
#include <givaro/udl.h>
#include "givaro/modular-extended.h"

#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/fflas/fflas_level1.inl"
#include "fflas-ffpack/fflas/fflas_level2.inl"
#include "fflas-ffpack/fflas/fflas_level3.inl"
#include "fflas-ffpack/fflas/fflas_enum.h"

namespace FFPACK {

    template<typename RNS>
    class RNSIntegerMod;
}
#include "fflas-ffpack/fflas/fflas_fscal_mp.inl"

#if defined(BENCH_PERF_FGEMM_MP) || defined(BENCH_PERF_TRSM_MP) || defined(BENCH_PERF_LQUP_MP)
#define BENCH_PERF_SCAL_MP
#define BENCH_MODP
#endif

namespace FFPACK {


    template<typename RNS>
    class RNSIntegerMod {
    public:
        typedef typename RNS::Element                   Element;
        typedef typename RNS::Element_ptr           Element_ptr;
        typedef typename RNS::ConstElement_ptr ConstElement_ptr;
        //typedef rnsRandIter<RNS> RandIter;
        class RandIter : public rnsRandIter<RNS> {
        private:
            const RNSIntegerMod<RNS> & _F;
            uint64_t _seed;
        public:
            RandIter(const RNSIntegerMod<RNS> &F, uint64_t seed=0) :rnsRandIter<RNS>(*F._rns,seed), _F(F), _seed(seed) {}

            typename RNS::Element& random(typename RNS::Element& elt) const {
                integer::seeding(_seed);
                integer tmp;
                integer::random_exact<true>(tmp,_F._p);
                _F.init(elt, tmp);
                return elt;
            }
        };



    protected:
        typedef typename RNS::BasisElement BasisElement;
        typedef Givaro::Modular<BasisElement> ModField;
        typedef Givaro::Integer integer;

        integer                              _p;
        std::vector<BasisElement, AlignedAllocator<BasisElement, Alignment::CACHE_LINE>>       _Mi_modp_rns;
        std::vector<BasisElement, AlignedAllocator<BasisElement, Alignment::CACHE_LINE>>       _iM_modp_rns;
        const RNS                         *_rns;
        Givaro::Modular<Givaro::Integer>     _F;
        RNSInteger<RNS>                  _RNSdelayed;
    public:
        Element                one, mOne,zero;

#ifdef BENCH_MODP
        mutable double t_modp, t_igemm, t_scal,t_trsm;
        mutable size_t n_modp;
#endif


        RNSIntegerMod(const integer& p, const RNS& myrns) : _p(p),
        _Mi_modp_rns(myrns._size*myrns._size),
        _iM_modp_rns((myrns._size+1)*myrns._size),
        _rns(&myrns),
        _F(p),
        _RNSdelayed(myrns){
            init(one,1);
            init(zero,0);
            init(mOne,-1);
            integer iM=0;
            size_t mysize=myrns._size;
            size_t mysizep1=myrns._size+1;
            integer sum=0;
            //std::cout << "M: " << myrns._M << std::endl;
            for (size_t i=0;i<mysize;i++){
                integer Mi = myrns._Mi[i] % _p;
                for (size_t j=0;j<mysize;j++){
                    _iM_modp_rns[i+j*mysizep1]=  iM % myrns._basis[j];
                    _Mi_modp_rns[i+j*mysize]=  Mi % myrns._basis[j];
                }
                iM+=myrns._M;iM%=_p;
                sum+=myrns._basis[i];
            }
            // last line of _iM_modp_rns corresponds to a quotient of _size
            for (size_t j=0;j<mysize;j++)
                _iM_modp_rns[mysize+j*mysizep1]=  iM % myrns._basis[j];
#ifdef BENCH_MODP
            t_modp=t_igemm=t_scal=t_trsm=0.;
            n_modp=0;
#endif
        }

        const rns_double& rns() const {return *_rns;}
        const RNSInteger<RNS>& delayed() const {return _RNSdelayed;}

        size_t size() const {return _rns->_size;}

        bool isOne(const Element& x) const {
            bool isone=true;
            for (size_t i=0;i<_rns->_size;i++)
                isone&= (one._ptr[i]== x._ptr[i]);
            return isone;
        }

        bool isMOne(const Element& x) const {
            bool ismone=true;
            for (size_t i=0;i<_rns->_size;i++)
                ismone&= (mOne._ptr[i]== x._ptr[i]);
            return ismone;
        }

        bool isZero(const Element& x) const {
            //write(std::cout,x)<<" == ";
            //write(std::cout,zero)<<std::endl;
            integer t1;
            t1=convert(t1,x)%_p;
            //std::cout<<"t1="<<t1<<std::endl;
            bool iszero=true;
            for (size_t i=0;i<_rns->_size;i++)
                iszero&= (zero._ptr[i]==x._ptr[i]);
            //std::cout<<(iszero || (t1==integer(0))?"zero":"nonzero")<<std::endl;
            return iszero || (t1==integer(0));
        }

        integer& characteristic(integer &p) const { return p=_p;}
        integer characteristic() const { return _p;}

        integer& cardinality(integer &p) const { return p=_p;}
        integer cardinality() const { return _p;}

            // just to test <0 in fflas_io
        integer minElement() const { return integer(0); }
            // just for bitsize
        integer maxElement() const { return _p-1; }

        Element& init(Element& x) const{
            if (x._ptr == NULL){
                x._ptr = FFLAS::fflas_new<BasisElement>(_rns->_size);
                x._stride=1;
                x._alloc=true;
            }
            return x;
        }
        Element& init(Element& x, const Givaro::Integer& y) const{
            init(x);
            size_t k =(_p.bitsize())/16+((_p.bitsize())%16?1:0);
            _rns->init(1,1,x._ptr,x._stride, &y,1,k);
            return x;
        }

        // assume this is the mod p operation
        Element& reduce (Element& x, const Element& y) const{
            Givaro::Integer tmp;
            convert(tmp,y);
            tmp %= _p;
            init (x,tmp);
            return x;
        }

        Element& reduce (Element& x) const{
            Givaro::Integer tmp;
            convert (tmp, x);
            tmp %= _p;
            return init (x, tmp);
        }

        Element& init(Element& x, const Element& y) const{
            return reduce (x, y);
        }


        Givaro::Integer convert(Givaro::Integer& x, const Element& y)const {
            _rns->convert(1,1,integer(0),&x,1,y._ptr,y._stride);
            return x;
        }

        Element& assign(Element& x, const Element& y) const {
            for(size_t i=0;i<_rns->_size;i++)
                x._ptr[i*x._stride] = y._ptr[i*y._stride];
            return x;
        }

        Element& add(Element& x, const Element& y, const Element& z) const {
            for(size_t i=0;i<_rns->_size;i++)
                _rns->_field_rns[i].add((x._ptr)[i*x._stride],
                                        (y._ptr)[i*y._stride],
                                        (z._ptr)[i*z._stride]);
            return x;
        }

        Element& sub(Element& x, const Element& y, const Element& z) const {
            for(size_t i=0;i<_rns->_size;i++)
                _rns->_field_rns[i].sub((x._ptr)[i*x._stride],
                                        (y._ptr)[i*y._stride],
                                        (z._ptr)[i*z._stride]);
            return x;
        }

        Element& neg(Element& x, const Element& y) const {
            for(size_t i=0;i<_rns->_size;i++)
                _rns->_field_rns[i].neg((x._ptr)[i*x._stride],
                                        (y._ptr)[i*y._stride]);
            return x;
        }

        Element& mul(Element& x, const Element& y, const Element& z) const {
            for(size_t i=0;i<_rns->_size;i++)
                _rns->_field_rns[i].mul((x._ptr)[i*x._stride],
                                        (y._ptr)[i*y._stride],
                                        (z._ptr)[i*z._stride]);
            return x;
        }


        Element& axpyin(Element& x, const Element& y, const Element& z) const {
            for(size_t i=0;i<_rns->_size;i++)
                _rns->_field_rns[i].axpyin((x._ptr)[i*x._stride],
                                           (y._ptr)[i*y._stride],
                                           (z._ptr)[i*z._stride]);
            return x;
        }

        Element& inv(Element& x, const Element& y) const {
            Givaro::Integer tmp;
            convert(tmp,y);
            _F.invin(tmp);
            init(x,tmp);
            return x;
        }

        bool areEqual(const Element& x, const Element& y) const {
            for(size_t i=0;i<_rns->_size;i++)
                if (!_rns->_field_rns[i].areEqual((x._ptr)[i*x._stride],(y._ptr)[i*y._stride]))
                    return false;
            return true;
        }
        std::ostream& write(std::ostream& os, const Element& y) const {
            // integer x;
            // convert(x,y);
            os<<"["<<(long) (y._ptr)[0];
            for(size_t i=1;i<_rns->_size;i++)
                os<<" , "<< (long) ((y._ptr)[i*y._stride]);
            return os<<" ]";
        }


        std::ostream& write(std::ostream& os) const {
            os<<" RNSIntegerMod("<<_p<<") with M:=[ "<< (long) _rns->_basis[0];
            for(size_t i=1;i<_rns->_size;i++)
                os<<" , "<< (long) _rns->_basis[i];
            return os<<" ]"<<std::endl;
        }


        void reduce_modp(size_t n, Element_ptr B) const{
#ifdef BENCH_MODP
            FFLAS::Timer chrono; chrono.start();
#endif

            size_t _size= _rns->_size;
            BasisElement *Gamma, *alpha, *A;
            A=B._ptr;
            size_t rda = B._stride;
            Givaro::ZRing<BasisElement> D;
            Gamma = FFLAS::fflas_new(D,_size,n);
            alpha = FFLAS::fflas_new(D,n);
            // compute Gamma
            typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
            FFLAS::fscal(_RNSdelayed, n, mmi, B, 1, typename RNS::Element_ptr(Gamma,n), 1);

            // compute A = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, _size, n, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, n, D.zero, A, rda);
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)_size, (int)n, (int)_size, 1.0 , _Mi_modp_rns.data(), (int)_size, Gamma, (int)n, 0, A, (int)rda);
#endif
            // compute alpha = _invbase.Gamma
            FFLAS::fgemv(D,FFLAS::FflasTrans, _size, n, D.one, Gamma, n, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);

            // compute ((z-(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
            for(size_t i=0;i<_size;i++){
                for(size_t j=0;j<n;j++){
                    size_t aa= (size_t)rint(alpha[j]);
#ifdef __FFLASFFPACK_DEBUG
                    if (aa>_size) {std::cout<<"RNS modp ERROR"<<std::endl;exit(1);}
#endif
                    A[j+i*rda]-=_iM_modp_rns[aa+i*(_size+1)];

                }
            }
            // reduce each row of A modulo m_i
            for (size_t i=0;i<_size;i++)
                FFLAS::freduce (_rns->_field_rns[i], n, A+i*rda, 1);

            FFLAS::fflas_delete(Gamma);
            FFLAS::fflas_delete(alpha);

#ifdef BENCH_MODP
            chrono.stop();
            t_modp+=chrono.usertime();
#endif
        }

        std::ostream& write_matrix(std::ostream& c,
                                   const double* E,
                                   int n, int m, int lda) const
        {
            c<<std::endl<<"***********************"<<std::endl;
            for (int i = 0; i<n;++i){
                for (int j=0; j<m;++j)
                    c << *(E+j+lda*i) << " ";
                c << std::endl;
            }
            c<<"***********************"<<std::endl;
            return c << std::endl;
        }
        std::ostream& write_matrix_long(std::ostream& c,
                                        const double * E,
                                        int n, int m, int lda) const
        {
            c<<std::endl<<"***********************"<<std::endl;
            for (int i = 0; i<n;++i){
                for (int j=0; j<m;++j)
                    c << (long)*(E+j+lda*i) << " ";
                c << std::endl;
            }
            c<<"***********************"<<std::endl;
            return c << std::endl;
        }

        void reduce_modp(size_t m, size_t n, Element_ptr B, size_t lda) const{
            const size_t mn=m*n;
            if (mn) {
#ifdef BENCH_MODP
            FFLAS::Timer chrono; chrono.start();
#endif
            //cout<<"REDUCE MOD WITH LDA!=N"<<endl;
            size_t _size= _rns->_size;
            BasisElement *Gamma, *alpha, *z, *A;
            A=B._ptr;
            size_t rda=B._stride;
            Gamma = FFLAS::fflas_new<BasisElement>(mn*_size);
            alpha = FFLAS::fflas_new<BasisElement>(mn);
            z     = FFLAS::fflas_new<BasisElement>(mn*_size);

            // compute Gamma
            typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
            FFLAS::fscal(_RNSdelayed, m, n, mmi, B, lda, typename RNS::Element_ptr(Gamma,mn), n);

            // compute Gamma = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
            Givaro::ZRing<BasisElement> D;
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,_size, mn, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, mn, D.zero, z, mn);
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)_size, (int)mn, (int)_size, 1.0 , _Mi_modp_rns.data(), (int)_size, Gamma, (int)mn, 0, z, (int)mn);
#endif
            // compute alpha = _invbase.Gamma
            FFLAS::fgemv(D, FFLAS::FflasTrans, _size, mn, D.one, Gamma, mn, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);

            // compute A=((Gamma--(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
            for(size_t k=0;k<_size;k++){
                for(size_t i=0;i<m;i++)
                    for(size_t j=0;j<n;j++){
                        size_t aa=(size_t)floor(alpha[j+i*n]+0.5);
#ifdef __FFLASFFPACK_DEBUG
                        if (aa>_size) {std::cout<<"RNS modp ERROR"<<std::endl;exit(1);}
#endif
                        A[j+i*lda+k*rda]= z[j+i*n+k*mn]-_iM_modp_rns[aa+k*(_size+1)];
                    }
            }

            // reduce each row of A modulo m_i
            for (size_t i=0;i<_size;i++)
                FFLAS::freduce (_rns->_field_rns[i], m, n, A+i*rda, lda);
            FFLAS::fflas_delete(Gamma);
            FFLAS::fflas_delete(alpha);
            FFLAS::fflas_delete(z);
#ifdef BENCH_MODP
            chrono.stop();
            t_modp+=chrono.usertime();
#endif
            }
        }

#ifdef __DLP_CHALLENGE

#define DELTA 27
        template<class SimdT>
        inline void splitSimd(const SimdT x, SimdT & x_h, SimdT & x_l) const {
            using simd = Simd<double>;
            using vect_t = typename simd::vect_t;
            vect_t vc = simd::set1((double)((1_ui64 << DELTA) + 1_ui64));
            vect_t tmp = simd::mul(vc, x);
            x_h = simd::add(tmp, simd::sub(x, tmp));
            x_l = simd::sub(x, x_h);
        }

        template<class SimdT>
        inline void multSimd(const SimdT va, const SimdT vb, SimdT & vs, SimdT & vt) const{
            using simd = Simd<double>;
            using vect_t = typename simd::vect_t;
            vect_t vah, val, vbh, vbl;
            splitSimd(va, vah, val);
            splitSimd(vb, vbh, vbl);
            vs = simd::mul(va, vb);
            vt = simd::add(simd::add(simd::sub(simd::mul(vah, vbh), vs), simd::mul(vah, vbl)), simd::add(simd::mul(val, vbh), simd::mul(val, vbl)));
        }

        template<class SimdT>
        inline SimdT multModSimd(const SimdT a, const SimdT b, const SimdT p, const SimdT ip, const SimdT np) const{
            using simd = Simd<double>;
            using vect_t = typename simd::vect_t;
            vect_t abh, abl, pqh, pql;
            multSimd(a, b, abh, abl);
            vect_t q = simd::floor(simd::mul(abh, ip));
            multSimd(p, q, pqh, pql);
            vect_t r = simd::add(simd::sub(abh, pqh), simd::sub(abl, pql));
            abh = simd::greater_eq(r, p);
            abl = simd::lesser(r, simd::zero());
            abh = simd::vand(abh, np);
            abl = simd::vand(abl, p);
            abh = simd::vor(abh, abl);
            return r = simd::add(r, abh);
        }

        inline void split(const double x, const int delta, double &x_h, double &x_l) const {
            double c = (double)((1_ui64 << delta) + 1_ui64);
            x_h = (c*x)+(x-(c*x));
            x_l = x - x_h;
        }

        inline void mult(const double a, const double b, double &s, double &t) const{
            double ah, al, bh, bl;
            s = a*b;
#ifdef __FMA__
            t = std::fma(-a, b, s);
#else
            split(a, DELTA, ah, al);
            split(b, DELTA, bh, bl);
            t = ((((-s+ah*bh)+(ah*bl))+(al*bh))+(al*bl));
#endif
        }

        inline double multmod(const double a, const double b, const double p, const double ip, const double np) const{
            double abh, abl, pqh, pql;
            mult(a, b, abh, abl);
            double q = floor(abh*ip);
            mult(p, q, pqh, pql);
            double r = (abh-pqh)+(abl-pql);
            if(r > p)
                r -= p;
            else if(r < 0)
                r += p;
            return r;
        }

        void reduce_modp_rnsmajor_scal_quad(size_t n, Element_ptr B) const {
            // std::cout << "modp scalar quad" << std::endl;
            // using namespace modp_details;
            using simd = Simd<BasisElement>;
            using vect_t = typename simd::vect_t;

            FFLAS::Timer T;
            size_t _size= _rns->_size;

            Givaro::ZRing<BasisElement> D;
            std::vector<Givaro::ModularExtended<double>> Fields;
            for(size_t i = 0 ; i < _size ; ++i){
                Fields.emplace_back(_rns->_basis[i]);
            }
            /*
               if((int64_t)B._ptr%simd::alignment == 0 && _size%simd::vect_size==0){
               for(size_t j = 0 ; j < n ; ++j){
               BasisElement *A, *Gamma, *tabTmp;
               A = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);
               Gamma = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);
               tabTmp = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);

               vect_t vA, vB, vG, vp, vnp, vip, vRNS, vtmp;

            // Compute Gamma
            for(size_t i = 0 ; i < _size ; i+= simd::vect_size){
            vB = simd::load(B._ptr+j*_size+i);
            vp = simd::load(_rns->_basis.data()+i);
            vip = simd::load(_rns->_invbasis.data()+i);
            vnp = simd::load(_rns->_negbasis.data()+i);
            vRNS = simd::load(_rns->_MMi.data()+i);
            vG = multModSimd(vB, vRNS, vp, vip, vnp);
            simd::store(Gamma+i, vG);
            }

            // Compute A=Gamma*Mi in rns
            for(size_t k = 0 ; k < _size ; ++k){
            for(size_t i = 0 ; i < _size ; i+= simd::vect_size){
            vG = simd::load(Gamma+i);
            vp = simd::set1(_rns->_basis[k]);
            vip = simd::set1(_rns->_invbasis[k]);
            vnp = simd::set1(_rns->_negbasis[k]);
            vRNS = simd::load(_Mi_modp_rns.data()+k*_size+i);
            vtmp = multModSimd(vG, vRNS, vp, vip, vnp);
            simd::store(tabTmp+i, vtmp);
            }
            for(size_t i = 0 ; i < _size ; ++i){
            Fields[k].addin(A[k], tabTmp[i]);
            }
            }
            double alpha = 0;
            for(size_t k = 0 ; k < _size ; ++k){
            alpha += Gamma[k]*_rns->_invbasis[k];
            }
            // -= alpha
            long aa= (long)rint(alpha);
            for(size_t k = 0; k < _size ; k++){
            Fields[k].sub(B._ptr[j*_size+k], A[k], _iM_modp_rns[aa+k*_size]);
            }
            FFLAS::fflas_delete(Gamma);
            FFLAS::fflas_delete(A);
            FFLAS::fflas_delete(tabTmp);
            }
            }else{
            for(size_t j = 0 ; j < n ; ++j){
            BasisElement *A, *Gamma, *tabTmp;
            A = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);
            Gamma = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);
            tabTmp = FFLAS::fflas_new<BasisElement>(_size, Alignment::CACHE_LINE);

            vect_t vA, vB, vG, vp, vnp, vip, vRNS, vtmp;

            // Compute Gamma
            size_t i = 0;
            for(; i < ROUND_DOWN(_size, simd::vect_size) ; i+= simd::vect_size){
            vB = simd::load(B._ptr+j*_size+i);
            vp = simd::load(_rns->_basis.data()+i);
            vip = simd::load(_rns->_invbasis.data()+i);
            vnp = simd::load(_rns->_negbasis.data()+i);
            vRNS = simd::load(_rns->_MMi.data()+i);
            vG = multModSimd(vB, vRNS, vp, vip, vnp);
            simd::store(Gamma+i, vG);
            }
            for(; i < _size ; ++i){
            Fields[i].mul(Gamma[i], B._ptr[j*_size+i], _rns->_MMi[i]);
    }

    // Compute A=Gamma*Mi in rns
    for(size_t k = 0 ; k < _size ; ++k){
        i = 0;
        A[k] = 0;
        for( ; i < ROUND_DOWN(_size, simd::vect_size); i+= simd::vect_size){
            vG = simd::load(Gamma+i);
            vp = simd::set1(_rns->_basis[k]);
            vip = simd::set1(_rns->_invbasis[k]);
            vnp = simd::set1(_rns->_negbasis[k]);
            vRNS = simd::load(_Mi_modp_rns.data()+k*_size+i);
            vtmp = multModSimd(vG, vRNS, vp, vip, vnp);
            simd::store(tabTmp+i, vtmp);
        }
        for(; i < _size ; ++i){
            Fields[k].mul(tabTmp[i], Gamma[i], _Mi_modp_rns[i]);
        }
        for(size_t i = 0 ; i < _size ; ++i){
            Fields[k].addin(A[k], tabTmp[i]);
        }
    }
    double alpha = 0;
    for(size_t k = 0 ; k < _size ; ++k){
        alpha += Gamma[k]*_rns->_invbasis[k];
    }
    // -= alpha
    long aa= (long)rint(alpha);
    for(size_t k = 0; k < _size ; k++){
        Fields[k].sub(B._ptr[j*_size+k], A[k], _iM_modp_rns[aa+k*_size]);
    }
    FFLAS::fflas_delete(Gamma);
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(tabTmp);
    }
        }
        //*/
        //*
#pragma omp parallel for schedule(static, 256)
        for(size_t i = 0 ; i < n; ++i){
            double* Ad;
            BasisElement *Gamma;
            Gamma = FFLAS::fflas_new<BasisElement>(_size);
            Ad = FFLAS::fflas_new<double>(_size);
            // Compute Gamma
            // std::cout << "B: " << std::endl;
            //  for(size_t j = 0 ; j < _size ; ++j){
            //  	std::cout << B._ptr[i*_size+j] << " ";
            //  }
            //  std::cout << std::endl;
            for(size_t k = 0; k < _size ; ++k){
                Fields[k].mul(Gamma[k], B._ptr[i*_size+k], _rns->_MMi[k]);
            }
            // std::cout << "Gamma: " << std::endl;
            // for(size_t j = 0 ; j < _size ; ++j){
            // 	std::cout << Gamma[j] << " ";
            // }
            // std::cout << std::endl;

            // std::cout << "MMi: " << std::endl;
            // for(size_t j = 0 ; j < _size ; ++j){
            // 	std::cout << _rns->_MMi[j] << " ";
            // }
            // std::cout << std::endl;

            // FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasTrans, n, _size, _size, D.one, Gamma, _size, _Mi_modp_rns.data(), _size, D.zero, A, _size);
            // Mul by Mi_modp
            for(size_t k = 0 ; k < _size ; ++k){
                Ad[k] = FFLAS::fdot(Fields[k], _size, Gamma, 1, _Mi_modp_rns.data()+k*_size,1);
            }
            // std::cout << "_Mi_modp_rns: " << std::endl;
            // std::cout << "[";
            // for(size_t j = 0 ; j < _size ; ++j){
            // 	std::cout << "[";
            // 	for(size_t k = 0 ; k < _size-1 ; ++k){
            // 		std::cout << _Mi_modp_rns[j*_size+k] << " ,";
            // 	}
            // 	std::cout << _Mi_modp_rns[j*_size+_size-1] << "],";
            // }
            // std::cout << "]" << std::endl;
            // std::cout << "Ad: " << std::endl;
            // for(size_t j = 0 ; j < _size ; ++j){
            // 	std::cout << Ad[j] << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "_iM_modp_rns: " << std::endl;
            // std::cout << "[";
            // for(size_t j = 0 ; j < _size ; ++j){
            // 	std::cout << "[";
            // 	for(size_t k = 0 ; k < _size-1 ; ++k){
            // 		std::cout << _iM_modp_rns[j*_size+k] << " ,";
            // 	}
            // 	std::cout << _iM_modp_rns[j*_size+_size-1] << "],";
            // }
            // std::cout << std::endl;

            // compute alpha
            // FFLAS::fgemv(D,FFLAS::FflasNoTrans, n, _size, D.one, Gamma, _size, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
            double alpha = 0;
            for(size_t k = 0 ; k < _size ; ++k){
                alpha += Gamma[k]*_rns->_invbasis[k];
            }
            // std::cout << "alpha: " << alpha << std::endl;
            // -= alpha
            // long aa= (int64_t)alpha;
            long aa= (long)rint(alpha);
            //std::cout << "aa: " << aa << std::endl;
            for(size_t k = 0; k < _size ; k++){
                // std::cout << Ad[k] << " - " << _iM_modp_rns[aa+k*_size] << " = ";
                Fields[k].sub(B._ptr[i*_size+k], Ad[k], _iM_modp_rns[aa+k*_size]);
                // std::cout<<B._ptr[i*_size+k]<<std::endl;
            }
            FFLAS::fflas_delete(Gamma);

            FFLAS::fflas_delete(Ad);
            // std::cout << std::endl;
            // std::cout << "====================================" << std::endl;
        }
        //*/
        // std::cout << std::endl;
        // _rns->reduce(n,B._ptr,1,true);
        }

#endif // __DLP_CHALLENGE

        void reduce_modp_rnsmajor(size_t n, Element_ptr B) const{
            // std::cout << "modp BLAS" << std::endl;
#ifdef BENCH_MODP
            FFLAS::Timer chrono; chrono.start();
#endif
            size_t _size= _rns->_size;
            BasisElement *Gamma, *alpha, *A;
            A=B._ptr;
            Givaro::ZRing<BasisElement> D;
            FFLAS::Timer T;
            // T.start();
            Gamma = FFLAS::fflas_new(D,n,_size);
            alpha = FFLAS::fflas_new(D,n,1);
            // T.stop();
            // std::cout << "Alloc: " << T << std::endl;
            // compute Gamma (NOT EFFICIENT)
            //for(size_t i=0;i<_size;i++)
            //
            // FFLAS::fscal(_rns->_field_rns[i], n, _rns->_MMi[i], A+i, _size, Gamma+i,_size);
            T.start();
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
            using simd = Simd<BasisElement>;
            using vect_t = typename simd::vect_t;

            if(((int64_t)A%simd::alignment == 0) && (_size%simd::vect_size==0)){
                auto MMi = _rns->_MMi.data();
                for(size_t i = 0 ; i < n ; ++i){
                    vect_t vA1, vA2, vMi1, vMi2, tmp1, tmp2, tmp3, v, max, basis, inv_, neg_;
                    size_t k = 0;
                    for( ; k < ROUND_DOWN(_size, simd::vect_size) ; k+=simd::vect_size){
                        basis = simd::load(_rns->_basis.data()+k);
                        inv_  = simd::load(_rns->_invbasis.data()+k);
                        max   = simd::load(_rns->_basisMax.data()+k);
                        neg_  = simd::load(_rns->_negbasis.data()+k);
                        vA1   = simd::load(A+i*_size+k);
                        vMi1  = simd::load(MMi+k);
                        v     = simd::mul(vA1, vMi1);
                        tmp1  = simd::floor(simd::mul(v, inv_));
                        tmp2  = simd::fnmadd(v, tmp1, basis);
                        tmp1  = simd::greater(tmp2, max);
                        tmp3  = simd::lesser(tmp2, simd::zero());
                        tmp1  = simd::vand(tmp1, neg_);
                        tmp3  = simd::vand(tmp3, basis);
                        tmp1  = simd::vor(tmp1, tmp3);
                        tmp2  = simd::add(tmp2, tmp1);
                        simd::store(Gamma+i*_size+k, tmp2);
                    }
                }
            }else{
                vect_t vA1, vA2, vMi1, vMi2, tmp1, tmp2, tmp3, v, max, basis, inv_, neg_;
                auto MMi = _rns->_MMi.data();
                for(size_t i = 0 ; i < n ; ++i){
                    size_t k = 0;
                    for( ; k < ROUND_DOWN(_size, simd::vect_size) ; k+=simd::vect_size){
                        basis = simd::load(_rns->_basis.data()+k);
                        inv_  = simd::load(_rns->_invbasis.data()+k);
                        max   = simd::load(_rns->_basisMax.data()+k);
                        neg_  = simd::load(_rns->_negbasis.data()+k);
                        vA1 = simd::loadu(A+i*_size+k);
                        vMi1 = simd::loadu(MMi+k);
                        v = simd::mul(vA1, vMi1);
                        tmp1  = simd::floor(simd::mul(v, inv_));
                        tmp2  = simd::fnmadd(v, tmp1, basis);
                        tmp1  = simd::greater(tmp2, max);
                        tmp3  = simd::lesser(tmp2, simd::zero());
                        tmp1  = simd::vand(tmp1, neg_);
                        tmp3  = simd::vand(tmp3, basis);
                        tmp1  = simd::vor(tmp1, tmp3);
                        tmp2  = simd::add(tmp2, tmp1);
                        simd::storeu(Gamma+i*_size+k, tmp2);
                    }
                    for(; k < _size ; ++k){
                        Gamma[i*_size+k] = A[i*_size+k]*MMi[k];
                        Gamma[i*_size+k] -= std::floor(Gamma[i*_size+k]*_rns->_invbasis[k])*_rns->_basis[k];
                        if(Gamma[i*_size+k] >= _rns->_basis[k]){
                            Gamma[i*_size+k] -= _rns->_basis[k];
                        }else if(Gamma[i*_size+k] < 0){
                            Gamma[i*_size+k] += _rns->_basis[k];
                        }
                    }
                }
            }
            // _rns->reduce(n,Gamma,1,true);
#else
            typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
            FFLAS::fscal(_RNSdelayed, n, mmi, B, 1, typename RNS::Element_ptr(Gamma,1), 1);
#endif
            T.stop();
            // std::cout << "Gamma: " << T << std::endl;


            // compute A = Gamma._Mi_modp_rns^T (note must be reduced mod m_i, but this is postpone to the end)
            T.start();
            FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasTrans, n, _size, _size, D.one, Gamma, _size, _Mi_modp_rns.data(), _size, D.zero, A, _size);
            T.stop();
            // std::cout << "fgemm: " << T << std::endl;
            // std::cout<<"fgemv (Y)...";
            //std::cout<<"fgemv (Y)..."<<n<<" -> "<<_size<<endl;;
            // compute alpha = Gamma._invbasis
            T.start();
            FFLAS::fgemv(D,FFLAS::FflasNoTrans, n, _size, D.one, Gamma, _size, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
            T.stop();
            // std::cout << "fgemv: " << T << std::endl;
            //std::cout<<"done"<<std::endl;
            T.start();
            // compute ((z-(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
            for(size_t j=0;j<n;j++){
                long aa= (long)rint(alpha[j]);
                for(size_t i=0;i<_size;i++){
                    //long aa=floor(alpha[j]+0.5);
                    A[j*_size+i]-=_iM_modp_rns[aa+i*_size];
                }
            }
            // vect_t viM;
            // auto iM = _iM_modp_rns.data();
            // for(size_t j = 0 ; j < n ; ++j){
            // 	long aa= (long)rint(alpha[j]);
            // 	for(int i = 0 ; i < ROUND_DOWN(_size, simd::vect_size) ; i+=simd::vect_size){
            // 		vA = simd::load(A+j*_size+i);
            // 		viM = simd::load(iM+aa)
            // 	}
            // }
            T.stop();
            // std::cout << "last: " << T << std::endl;
            T.start();
            // reduce each column of A modulo m_i
            _rns->reduce(n,A,1,true);
            T.stop();
            // std::cout << "reduce: "<< T << std::endl;

            // T.start();
            FFLAS::fflas_delete(Gamma);
            FFLAS::fflas_delete(alpha);
            // T.stop();
            // std::cout << "delete: " << T << std::endl;
#ifdef BENCH_MODP
            chrono.stop();
            t_modp+=chrono.usertime();
#endif

        }


    }; // end of class RNSIntegerMod

} // end of namespace FFPACK


namespace FFLAS {

    // specialization for the fflas alloc function
    template<>
    inline FFPACK::rns_double_elt_ptr
    fflas_new(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m, const Alignment align){
        return FFPACK::rns_double_elt_ptr(FFLAS::fflas_new<double>(m*F.size(),align),m);
    }

    template<>
    inline FFPACK::rns_double_elt_ptr
    fflas_new(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m,  const size_t n, const Alignment align){
        return fflas_new(F, m*n, align);
    }

    // function to convert from integer to RNS (note: this is not the finit function from FFLAS, extra k)
    template<typename RNS>
    void finit_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n, size_t k,
                   const Givaro::Integer *B, const size_t ldb, typename RNS::Element_ptr A)
    {
        F.rns().init(m,n,A._ptr,A._stride, B,ldb,k);
    }
    template<typename RNS>
    void finit_trans_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n, size_t k,
                         const Givaro::Integer *B, const size_t ldb, typename RNS::Element_ptr A)
    {
        F.rns().init_transpose(m,n,A._ptr,A._stride, B,ldb,k);
    }

    // function to convert from RNS to integer (note: this is not the fconvert function from FFLAS, extra alpha)
    template<typename RNS>
    void fconvert_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n,
                      Givaro::Integer alpha, Givaro::Integer *B, const size_t ldb, typename RNS::ConstElement_ptr A)
    {
        F.rns().convert(m,n,alpha,B,ldb,A._ptr,A._stride);
    }
    template<typename RNS>
    void fconvert_trans_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n,
                            Givaro::Integer alpha, Givaro::Integer *B, const size_t ldb, typename RNS::ConstElement_ptr A)
    {
        F.rns().convert_transpose(m,n,alpha,B,ldb,A._ptr,A._stride);
    }

} // end of namespace FFLAS
#undef DELTA
#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
