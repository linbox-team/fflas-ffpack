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
//                        Test for nullspace
//--------------------------------------------------------------------------
// Authors: Clément Pernet
//			Alexis Breust (clean-up)
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>

#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/timer.h"

template <class Field>
std::string checkingMessage(const Field& F)
{
    std::ostringstream message;
    std::ostringstream fieldNameStream;
    F.write(fieldNameStream);

    message << "Checking ";
    message.fill('.');
    message.width(40);
    message << fieldNameStream.str() << " ... ";

    return message.str();
}

/**
 * If file is not empty, read it and set m, n, lda and r.
 * Otherwise, generate a random matrix of size m x n with random lda.
 */
template <class Field>
typename Field::Element_ptr readOrRandomMatrixWithRankAndRandomRPM(const Field& F, std::string file, size_t& m, size_t& n,
                                                                   size_t& lda, size_t& r, uint64_t seed)
{
    typename Field::Element_ptr A = nullptr;

    if (!file.empty()) {
        FFLAS::ReadMatrix(file, F, m, n, A);
        lda = n;
        r = FFPACK::Rank(F, m, n, A, lda);
    }
    else {
        lda = std::max(m, n) + (rand() % 13);
        A = FFLAS::fflas_new(F, m, lda);
        typename Field::RandIter G(F, seed);
        FFPACK::RandomMatrixWithRankandRandomRPM(F, m, n, r, A, lda, G);
    }

    return A;
}

template <class Field>
bool test_nullspace(Field& F, FFLAS::FFLAS_SIDE side, size_t m, size_t n, size_t r, typename Field::Element_ptr A, size_t lda)
{
    // @note As NullSpaceBasis mutates A, we make a copy of it.
    auto ACopy = FFLAS::fflas_new(F, m, lda);
    FFLAS::fassign(F, m, n, A, lda, ACopy, lda);

    size_t ldns = 0u;
    size_t NSdim = 0u;
    typename Field::Element_ptr NS;
    FFPACK::NullSpaceBasis(F, side, m, n, ACopy, lda, NS, ldns, NSdim);
    FFLAS::fflas_delete(ACopy);

#if defined(__FFLAS_FFPACK_DEBUG)
    std::cout << std::endl;
    std::cout << "A: " << m << "x" << n << " (rank " << r << ")" << std::endl;
    std::cout << "NS: " << NSdim << std::endl;
#endif
    if (side == FFLAS::FFLAS_SIDE::FflasRight) {
        // Right nullspace dimension + Rank == Matrix column dimension
        if (NSdim + r != n) return false;
        // Ensure nullspace is full rank
        auto NSCopy = FFLAS::fflas_new(F, n, NSdim);
        FFLAS::fassign(F, n, NSdim, NS, NSdim, NSCopy, NSdim);
        size_t rank = FFPACK::Rank(F, n, NSdim, NSCopy, NSdim);
        FFLAS::fflas_delete(NSCopy);
        if (rank != NSdim) return false;

        // Check that NS is a nullspace
        auto C = FFLAS::fflas_new(F, m, NSdim);
        FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, NSdim, n, F.one, A, lda, NS, ldns, F.zero, C, NSdim);
        if (!FFLAS::fiszero(F, m, NSdim, C, NSdim)) return false;
        FFLAS::fflas_delete(C);
    }
    else {
        // Left nullspace dimension + Rank == Matrix row dimension
        if (NSdim + r != m) return false;

        // Ensure nullspace is full rank
        auto NSCopy = FFLAS::fflas_new(F, NSdim, m);
        FFLAS::fassign(F, NSdim, m, NS, m, NSCopy, m);
        size_t rank = FFPACK::Rank(F, NSdim, m, NSCopy, m);
        FFLAS::fflas_delete(NSCopy);
        if (rank != NSdim) return false;

        // Check that NS is a nullspace
        auto C = FFLAS::fflas_new(F, NSdim, n);
        FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, NSdim, n, m, F.one, NS, ldns, A, lda, F.zero, C, n);
        if (!FFLAS::fiszero(F, NSdim, n, C, n)) return false;
        FFLAS::fflas_delete(C);
    }

    FFLAS::fflas_delete(NS);
    return true;
}

template <class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t r, size_t iters, std::string file, uint64_t& seed)
{
    bool ok = true;

    for (auto i = 0u; ok && i < iters; ++i) {
        Field* F = FFPACK::chooseField<Field>(q, b, seed);
        if (F == nullptr) return true;
        std::cout << checkingMessage(*F);

        size_t lda = 0u;
        auto A = readOrRandomMatrixWithRankAndRandomRPM(*F, file, m, n, lda, r, seed++);

        // The test indeed
        ok = ok && test_nullspace(*F, FFLAS::FFLAS_SIDE::FflasLeft, m, n, r, A, lda);
        ok = ok && test_nullspace(*F, FFLAS::FFLAS_SIDE::FflasRight, m, n, r, A, lda);

        FFLAS::fflas_delete(A);
        delete F;

        std::cout << (ok ? "PASSED" : "FAILED") << std::endl;
    }

    return ok;
}

int main(int argc, char** argv)
{
    Givaro::Integer q = -1;
    size_t b = 0;
    size_t m = 100;
    size_t n = 120;
    size_t r = 70;
    size_t iters = 3;
    uint64_t seed = FFLAS::getSeed();
    bool loop = false;
    std::string file;

    Argument as[] = {{'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER, &q},
        {'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT, &b},
        {'m', "-m M", "Set the row dimension of the matrix.", TYPE_INT, &m},
        {'n', "-n N", "Set the column dimension of the matrix.", TYPE_INT, &n},
        {'r', "-r R", "Set the rank.", TYPE_INT, &r},
        {'i', "-i I", "Set number of iterations.", TYPE_INT, &iters},
        {'s', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed},
        {'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop},
        {'f', "-f file", "Set input file", TYPE_STR, &file},
        END_OF_ARGUMENTS};

    FFLAS::parseArguments(argc, argv, as);

    if (r > std::min(m, n)) {
        r = std::min(m, n);
    }

    srand(seed);

    bool ok = true;
    do {
        auto lastKnownSeed = seed;
        ok = ok && run_with_field<Givaro::Modular<float>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::Modular<double>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::Modular<int32_t>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::Modular<int64_t>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<float>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<double>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int32_t>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int64_t>>(q, b, m, n, r, iters, file, seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer>>(q, 5, m / 6, n / 6, r / 6, iters, file, seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer>>(q, (b ? b : 512), m / 6, n / 6, r / 6, iters, file, seed);
        if (!ok) std::cerr << "Failed with seed: " << lastKnownSeed << std::endl;
    } while (loop && ok);

    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
