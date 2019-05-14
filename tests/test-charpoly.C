/*
 * Copyright (C) FFLAS-FFPACK
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

//--------------------------------------------------------------------------
//                        Test for charpoly
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

//#define ENABLE_ALL_CHECKINGS


#include <iostream>
#include <iomanip>
#include "givaro/modular.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/ffpack/ffpack.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

#include <random>
#include <chrono>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

template<class Field, class RandIter>
bool launch_test(const Field & F, size_t n, typename Field::Element * A, size_t lda,
                 size_t nbit, RandIter& G, FFPACK::FFPACK_CHARPOLY_TAG CT)
{
    std::ostringstream oss;
    switch (CT){
    case FfpackAuto: oss<<"Automated variant choice"; break;
    case FfpackLUK: oss<<"LUKrylov variant"; break;
    case FfpackKG: oss<<"Keller-Gehrig variant"; break;
    case FfpackDanilevski: oss<<"Danilevskii variant"; break;
    case FfpackKGFast:  oss<<"KGFast variant"; break;
    case FfpackKGFastG: oss<<"KGFastG variant"; break;
    case FfpackHybrid: oss<<"Hybrid variant"; break;
    case FfpackArithProgKrylovPrecond: oss<<"Precond. ArithProg variant"; break;
    default: oss<<"LUKrylov variant"; break;
    }
    F.write(oss<<" over ");
    std::cout.fill('.');
    std::cout<<"Checking ";
    std::cout.width(70);
    std::cout<<oss.str();
    std::cout<<"...";

    typedef typename Givaro::Poly1Dom<Field> PolRing;
    typedef typename PolRing::Element Polynomial;
    Polynomial charp(n+1);

    typename Field::Element_ptr B = FFLAS::fflas_new(F, n,n);
    FFLAS::fassign(F, n, n, A, lda, B, n);

    Checker_charpoly<Field,Polynomial> checker(F,n,A,lda);

    PolRing R(F);

    FFPACK::CharPoly (R, charp, n, A, lda, G, CT);

    try{
        checker.check(charp);
    }
    catch (FailureCharpolyCheck){
        std::cerr<<"CharPoly Checker FAILED"<<std::endl;
        FFLAS::fflas_delete (B);
        return false;
    }

    FFLAS::fassign(F, n, n, B, n, A, lda);

    // Checking trace(A) == charp[n-1]
    typename Field::Element trace;
    F.init(trace);
    F.assign(trace, F.zero);
    for (size_t i = 0; i < n; i++)
        F.subin (trace, B [i*(n+1)]);
    if (!F.areEqual(trace, charp[n-1])){
        std::cerr<<"FAILED: trace = "<<trace<<" P["<<n-1<<"] = "<<charp[n-1]<<std::endl;
        std::cerr<<" P = "<<charp<<std::endl;
        FFLAS::fflas_delete (B);
        return false;
    }

    // Checking det(A) == charp[0]
    typename Field::Element det;
    F.init(det);
    FFPACK::Det(F, det, n, B, n);
    FFLAS::fflas_delete (B);

    if (n&1) F.negin(det); // p0 == (-1)^n det

    if (!F.areEqual(det,charp[0])){
        std::cerr<<"FAILED: det = "<<det<<" P["<<0<<"] = "<<charp[0]<<std::endl;
        std::cerr<<" P = "<<charp<<std::endl;
        return false;
    }

    std::cout<<"PASSED"<<std::endl;
    return true ;
}

template<class Field>
bool run_with_field(const Givaro::Integer p, uint64_t bits, size_t n, std::string file, int variant, size_t iter, uint64_t seed){
    FFPACK::FFPACK_CHARPOLY_TAG CT;
    switch (variant){
    case 0: CT = FfpackAuto; break;
    case 1: CT = FfpackDanilevski; break;
    case 2: CT = FfpackLUK; break;
    case 3: CT = FfpackArithProgKrylovPrecond; break;
    case 4: CT = FfpackKG; break;
    case 5: CT = FfpackKGFast; break;
    case 6: CT = FfpackHybrid; break;
    case 7: CT = FfpackKGFastG; break;
    default: CT = FfpackAuto; break;
    }

    bool passed = true;
    size_t lda = n;

    for (size_t i=0;i<iter;i++){
        Field* F= chooseField<Field>(p,bits,seed);
        if (F==nullptr){
            return true;
        }
        typename Field::RandIter R(*F,seed++);

        typename Field::Element * A=NULL;

        if (!file.empty()) {
            /* user provided test matrix */
            //const char * filestring = file.c_str();
            FFLAS::ReadMatrix(file,*F,n,n,A);
        } else {
            /* Random matrix test */
            A = FFLAS::fflas_new(*F,n,lda);
            A = FFPACK::RandomMatrix(*F,n,n,A,lda,R);
        }

        if (variant) // User provided variant
            passed = passed && launch_test<Field>(*F, n, A, lda, iter, R, CT);
        else{ // No variant specified, testing them all
            passed = passed && launch_test<Field>(*F, n, A, lda, iter, R, FfpackDanilevski);
            passed = passed && launch_test<Field>(*F, n, A, lda, iter, R, FfpackLUK);
            passed = passed && launch_test<Field>(*F, n, A, lda, iter, R, FfpackArithProgKrylovPrecond);
            passed = passed && launch_test<Field>(*F, n, A, lda, iter, R, FfpackAuto);
            //passed = passed && launch_test<Field>(F, n, A, lda, iter, FfpackKG); // fails (variant only implemented for benchmarking
            //passed = passed && launch_test<Field>(*F, n, A, lda, iter, FfpackKGFast); // generic: does not work with any matrix
            //passed = passed && launch_test<Field>(*F, n, A, lda, iter, FfpackKGFastG); // generic: does not work with any matrix
            //passed = passed && launch_test<Field>(*F, n, A, lda, iter, FfpackHybrid); // fails with small characteristic
        }
        FFLAS::fflas_delete (A);
        delete F;
    }
    if (!passed) std::cerr<<std::endl<<"Failed with seed = "<<seed<<std::endl;
    return passed;
}

int main(int argc, char** argv)
{
    Givaro::Integer q = -1; // characteristic
    uint64_t     bits = 0;       // bit size
    size_t       iter = 2; // repetitions
    size_t       n = 150;
    std::string  file = "" ; // file where
    bool loop = false; // loop infintely
    std::string  mat_file = "" ; // input matrix file
    int variant = 0; // default value 0: test all variants
    uint64_t seed = getSeed();

    std::cout<<setprecision(17);

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic.", TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random elements.",         TYPE_INT , &bits},
        { 'n', "-n N", "Set the size of the matrix.", TYPE_INT , &n },
        { 'i', "-i I", "Set number of repetitions.", TYPE_INT , &iter },
        { 'f', "-f file", "Set input file", TYPE_STR, &mat_file },
        { 'a', "-a algorithm", "Set the algorithm variant", TYPE_INT, &variant },
        { 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);

    srand(seed);
    bool passed = true;
    do {
        passed = passed && run_with_field<Givaro::ModularBalanced<float> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<float> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::ModularBalanced<double> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<double> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::ModularBalanced<int32_t> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<int32_t> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::ModularBalanced<int64_t> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<int64_t> >(q, bits, n, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<Givaro::Integer> >(q, 6, n/2, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::Modular<Givaro::Integer> >(q, (bits?bits:512), n/4, mat_file, variant, iter, seed);
        passed = passed && run_with_field<Givaro::ZRing<Givaro::Integer> >(q, (bits?bits:80_ui64), n/4, mat_file, variant, iter, seed);

        // if ((i+1)*100 % nbit == 0)
        // 	std::cerr<<double(i+1)/nbit*100<<" % "<<std::endl;
    } while (loop && passed);
    return !passed;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
