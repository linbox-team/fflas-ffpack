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

#ifndef __FFPACK_ludivine_mp_INL
#define __FFPACK_ludivine_mp_INL

#include <givaro/modular-integer.h>
#include <givaro/givinteger.h>

#ifdef BENCH_PERF_LQUP_MP
#define BENCH_PERF_FGEMM_MP
#endif

#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack_ludivine.inl"

namespace FFPACK {

    template <>
    inline size_t
    LUdivine (const Givaro::Modular<Givaro::Integer>& F,
              const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
              const size_t M, const size_t N,
              typename Givaro::Integer* A, const size_t lda,
              size_t*P, size_t *Q,
              const FFPACK::FFPACK_LU_TAG LuTag,
              const size_t cutoff
             )
    {
        const size_t K = std::max(M,N);
        if (K) {
#ifdef BENCH_PERF_LQUP_MP
        double t_init=0, t_lqup=0, t_mod=0, t_rec=0;
        FFLAS::Timer chrono;
        chrono.start();
#endif
        Givaro::Integer p;
        F.cardinality(p);
        size_t logp=p.bitsize();

        // compute bit size of feasible prime
        size_t _k=std::max(K+1,logp/20), lk=0;
        while ( _k ) {_k>>=1; ++lk;}
        size_t prime_bitsize= (53-lk)>>1;

        // construct rns basis
        Givaro::Integer maxC= (p-1)*(p-1)*(p-1)*uint64_t(K);
        uint64_t n_pr =uint64_t(ceil(double(maxC.bitsize())/double(prime_bitsize)));
        maxC=(p-1)*(p-1)*uint64_t(K)*(1<<prime_bitsize)*n_pr;


        FFPACK::rns_double RNS(maxC, prime_bitsize, true);
        FFPACK::RNSIntegerMod<FFPACK::rns_double> Zp(p, RNS);
#ifdef BENCH_PERF_LQUP_MP
        chrono.stop();
        t_init+=chrono.usertime();
        chrono.clear();chrono.start();
#endif
        // compute A in RNS
        FFPACK::rns_double::Element_ptr Ap;
        Ap = FFLAS::fflas_new(Zp,M,N);
        FFLAS::finit_rns(Zp,M,N,(logp/16)+(logp%16?1:0),A,lda,Ap);


#ifdef BENCH_PERF_LQUP_MP
        chrono.stop();
        t_mod+=chrono.usertime();
        chrono.clear();chrono.start();
#endif

        // call lqup in rns
        size_t R=FFPACK::LUdivine(Zp, Diag, trans, M, N, Ap, N, P, Q, LuTag, cutoff);

        //std::cout<<"LUDivine RNS done"<<std::endl;
#ifdef BENCH_PERF_LQUP_MP
        chrono.stop();
        t_lqup+=chrono.usertime();
        chrono.clear();chrono.start();
#endif
        // reconstruct the result
        FFLAS::fconvert_rns(Zp,M,N,F.zero,A,lda,Ap);
#ifdef BENCH_PERF_LQUP_MP
        chrono.stop();
        t_rec+=chrono.usertime();
        chrono.clear();chrono.start();
#endif
        // reduce it modulo p
        FFLAS::freduce (F,M,N,A,lda);

#ifdef BENCH_PERF_LQUP_MP
        chrono.stop();
        //t_rec+=chrono.usertime();
        std::cout<<"LUDIVINE RNS PERF:"<<std::endl;
        std::cout<<"  ---  RNS basis size: "<<Zp.size() <<std::endl;
        std::cout<<"  ***      init  : "<<t_init<<std::endl;
        std::cout<<"  ***  rns  mod  : "<<t_mod<<std::endl;
        std::cout<<"  ***  rns lqup  : "<<t_lqup<<" ( igemm="<<Zp.t_igemm<<" ftrsm="<<Zp.t_trsm<<" scal="<<Zp.t_scal
        <<" modp="<<Zp.t_modp<<std::endl;
        std::cout<<"  ***  rns  rec  : "<<t_rec<<std::endl;
        std::cout<<"  ***       mod  : "<<chrono.usertime()<<std::endl;

#endif
        FFLAS::fflas_delete(Ap);
        return R;
        } else {
            return 0;
        }
    }

} // namespace FFPACK

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
