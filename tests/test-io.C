/*
 * Copyright (C) the FFLAS-FFPACK group 2017
 * Written by Clément Pernet
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


#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <random>

#include <givaro/modular.h>
#include <givaro/zring.h>
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
using namespace FFLAS;
using Givaro::Modular;

template<class Element> struct CompactElement{typedef Element type;};
template<> struct CompactElement<double> {typedef int32_t type;};
template<> struct CompactElement<float>  {typedef int16_t type;};
template<> struct CompactElement<int64_t>{typedef int32_t type;};
template<> struct CompactElement<int32_t>{typedef int16_t type;};
template<> struct CompactElement<int16_t>{typedef int8_t type;};

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t iters, uint64_t seed){

    bool ok=true;
    int nbit =(int) iters;
    while (ok && nbit){
        Field* F = FFPACK::chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
        typename Field::RandIter G(*F, seed++);
        std::ostringstream oss;
        F->write(oss);

        std::cout.fill('.');
        std::cout<<"Checking ";


        string file_dense = "data/mat.dense";
        string file_sms = "data/mat.sms";
        string file_binary = "data/mat.bin";
        string file_compact_binary = "data/mat.cbin";

        typename Field::Element_ptr A = fflas_new (*F, m, n);
        FFPACK::RandomMatrix (*F, m, n, A, n, G);
        typename Field::Element_ptr B;

        // Testing Dense format
        WriteMatrix (file_dense,*F,m,n,A,n, FflasDense);
        ReadMatrix (file_dense,*F,m,n,B, FflasDense);
        ok = ok && fequal (*F, m, n, A, n, B, n);
        if (ok) oss<<" Dense (ok)";
        else oss<<" Dense (KO)"<<std::endl;
        fflas_delete(B);

        // Testing SMS format
        WriteMatrix (file_sms,*F,m,n,A,n, FflasSMS);
        ReadMatrix (file_sms,*F,m,n,B, FflasSMS);
        ok = ok && fequal (*F, m, n, A, n, B, n);
        if (ok) oss<<" SMS (ok)";
        else oss<<" SMS (KO)";
        fflas_delete(B);

        // Testing Binary format
        WriteMatrix (file_binary,*F,m,n,A,n, FflasBinary);
        ReadMatrix (file_binary,*F,m,n,B, FflasBinary);
        ok = ok && fequal (*F, m, n, A, n, B, n);
        if (ok) oss<<" Bin (ok)";
        else oss<<" Bin (KO)";
        fflas_delete(B);
        // Testing compact Binary format
        typedef Givaro::ZRing<typename CompactElement<typename Field::Element>::type> CompactField;
        CompactField Z;
        typename CompactField::Element_ptr Az = fflas_new(Z,m,n);
        B = fflas_new(*F,m,n);
        typename CompactField::Element_ptr  Bz = NULL;
        fconvert(*F,m,n,Az,n,A,n);
        WriteMatrix (file_compact_binary,Z,m,n,Az,n, FflasBinary);
        ReadMatrix (file_compact_binary,Z,m,n,Bz, FflasBinary);
        finit (*F,m,n,Bz,n,B,n);
        ok = ok && fequal (*F, m, n, A, n, B, n);
        if (ok) oss<<" Compact Bin (ok)";
        else oss<<" Compact Bin (KO)";
        fflas_delete(Az);
        fflas_delete(B);
        fflas_delete(Bz);

        // Testing Autodetection of Binary format
        ReadMatrix (file_binary,*F,m,n,B, FflasAuto);
        ok = ok && fequal (*F, m, n, A, n, B, n);
        if (ok) oss<<" Auto Bin (ok)";
        else oss<<" Auto Bin (KO)";
        fflas_delete(B);

        std::cout.width(75);
        std::cout<<oss.str();
        std::cout<<" ... ";

        if (ok) std::cout << "PASSED"<<std::endl;
        else std::cout << "FAILED"<<std::endl;
        fflas_delete(A);
        delete F;
        nbit--;
    }
    return 0;
}


int main(int argc, char** argv){
    cerr<<setprecision(20);
    Givaro::Integer q=-1;
    size_t b=0;
    size_t m=53;
    size_t n=97;
    size_t iters=3;
    bool loop=false;
    uint64_t seed=getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'm', "-m M", "Set the row dimension of the matrix.",      TYPE_INT , &m },
        { 'n', "-n N", "Set the column dimension of the matrix.", TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
    srand(seed);
    bool ok=true;
    do{
        run_with_field<Modular<double> >(q,b,m,n,iters,seed);

    } while(loop && ok);
    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
