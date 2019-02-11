/*
 * Copyright (C) FFLAS-FFPACK
 * This file is Free Software and part of FFLAS-FFPACK.
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
#include <givaro/modular.h>
#include <recint/rint.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include <stdlib.h>
#include <stdio.h>

// This is required to compute min (max<size_t>, x) and return a size_t
#define MAX_WITH_SIZE_T(x) ( (static_cast<uint64_t>(std::numeric_limits<size_t>::max()) <  x)? std::numeric_limits<size_t>::max() : x )


template <class Field>
bool test (Givaro::Integer p, size_t kmax){
    Field F(p);
    FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Auto, FFLAS::ModeCategories::DelayedTag> MMH(F, 0);
    size_t k = MMH.MaxDelayedDim(0);
    if (kmax!=k)
        F.write(std::cerr)<<": expected: "<<kmax<<" got: "<<k<<std::endl;
    return kmax == k;
}
int main() {

    bool ok=true;

    // kmax = floor(2^53 / (p-1)^2)
    ok = ok && test<Givaro::Modular<double>  >(17, MAX_WITH_SIZE_T(35184372088831ULL));
    ok = ok && test<Givaro::Modular<double>  >(65521, 2098176);
    ok = ok && test<Givaro::Modular<double>  >(67108859,2);
    // kmax = floor(2^53 / ((p-1)/2)^2)
    ok = ok && test<Givaro::ModularBalanced<double>  >(17,MAX_WITH_SIZE_T(140737488355327ULL));
    ok = ok && test<Givaro::ModularBalanced<double>  >(65521,8392705);
    ok = ok && test<Givaro::ModularBalanced<double>  >(67108859,8);
    // kmax = floor(2^24 / (p-1)^2)
    ok = ok && test<Givaro::Modular<float> > (17,65535);
    ok = ok && test<Givaro::Modular<float> > (2039,4);
    // kmax = floor(2^24 / ((p-1)/2)^2)
    ok = ok && test<Givaro::ModularBalanced<float> >(17,262143);
    ok = ok && test<Givaro::ModularBalanced<float> > (2039,16);

    // kmax = floor(2^53 / (p-1)^2)
    ok = ok && test<Givaro::Modular<int64_t>  >(17, MAX_WITH_SIZE_T(36028797018963967));
    ok = ok && test<Givaro::Modular<int64_t>  >(65521, MAX_WITH_SIZE_T(2148532608));
    ok = ok && test<Givaro::Modular<int64_t>  >(1147482977,7);
    // kmax = floor(2^53 / ((p-1)/2)^2)
    ok = ok && test<Givaro::ModularBalanced<int64_t>  >(17, MAX_WITH_SIZE_T(144115188075855871));
    ok = ok && test<Givaro::ModularBalanced<int64_t>  >(65521, MAX_WITH_SIZE_T(8594130432));
    ok = ok && test<Givaro::ModularBalanced<int64_t>  >(1147482977,28);

    // kmax = floor(2^31 / (p-1)^2)
    ok = ok && test<Givaro::Modular<int32_t>  >(17,8388607);
    ok = ok && test<Givaro::Modular<int32_t>  >(24571,3);
    // kmax = floor(2^31 / ((p-1)/2)^2)
    ok = ok && test<Givaro::ModularBalanced<int32_t>  >(17,33554431);
    ok = ok && test<Givaro::ModularBalanced<int32_t>  >(24571,14);

    // kmax = maxsize_t
    ok = ok && test<Givaro::Modular<Givaro::Integer>  >(17, std::numeric_limits<size_t>::max());
    ok = ok && test<Givaro::Modular<Givaro::Integer>  >(Givaro::Integer("46768052394588893382517914646921056628989841375373"),std::numeric_limits<size_t>::max());
    // kmax = maxsize_t
    ok = ok && test<Givaro::Modular<RecInt::rint<8> > >(17, std::numeric_limits<size_t>::max());
    ok = ok && test<Givaro::Modular<RecInt::rint<8> > >(Givaro::Integer("166153499473114484112975882535042793"),2097152);
    return !ok;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
