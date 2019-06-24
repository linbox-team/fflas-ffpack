/*
 * Copyright (C) FFLAS-FFPACK
 * Written by David Lucas
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
//                        Test for minpoly
//
//--------------------------------------------------------------------------
// David Lucas
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/test-utils.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/modular-integer.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givpoly1.h>

using namespace std;
using namespace FFLAS;
using namespace FFPACK;
using Givaro::Modular;
using Givaro::ModularBalanced;


template<typename Field, class RandIter>
bool check_minpoly(const Field &F, size_t n, RandIter& G)
{
    typedef typename Field::Element_ptr Element_ptr;
    typedef vector<typename Field::Element> Polynomial;
    size_t lda, ldv;
    Element_ptr A, V, Vcst;

    //Default
    lda = n;
    ldv = n;

    /*Create variables used for testing (matrices, vectors and polynomials) */

    A = FFLAS::fflas_new(F, n, lda);
    V = FFLAS::fflas_new(F, n+1, ldv);
    Vcst = FFLAS::fflas_new(F, n);
    Polynomial minP;


    FFPACK::RandomMatrix (F, n, n, A, lda, G);

    FFPACK::NonZeroRandomMatrix(F, 1, n, V, ldv, G);
    FFLAS::fassign(F, n, V, 1, Vcst, 1); //MatVecMinPoly modifies V, we store it in Vcst beforehand

    FFPACK::MatVecMinPoly(F, minP, n, A, lda, V, 1);
    FFLAS::fflas_delete(V);

    /*Check that minP is monic*/

    size_t deg = minP.size() - 1;
    if(!(F.areEqual(minP[deg], F.one)))
        return false;

    // Krylov matrix computation
    size_t ldk = n;
    Element_ptr K = FFLAS::fflas_new(F, deg+1, ldk);
    FFLAS::fassign(F, n, Vcst, 1, K, 1);
    FFLAS::fflas_delete(Vcst);
    Element_ptr Kptr = K;
    for(size_t i = 0; i < deg; ++i, Kptr += ldk)
        FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, lda, Kptr, 1, F.zero, Kptr+ldk, 1);


    /*Check that minP(A).V is zero*/
    Element_ptr E = FFLAS::fflas_new(F, n);
    FFLAS::fzero(F, n, E, 1);

    for(size_t i = 0; i < deg+1; ++i)
        FFLAS::faxpy(F, n, minP[i], K+i*ldk, 1, E, 1);

    if (!FFLAS::fiszero(F, n, E, 1))
    {
        cout<<"NONZEROERROR"<<endl;
        FFLAS::fflas_delete(E);
        return false;
    }

    FFLAS::fflas_delete(E);


    /* Check minimality of minP */



    // minP factorization
    typedef Givaro::Poly1FactorDom<Field, Givaro::Dense> PolyDom; //defines a polynomial domain for Givaro
    typedef typename PolyDom::Element FieldPoly; //defines an element over this polynomial domain (casting purposes)
    vector<FieldPoly> factors;
    vector<uint64_t> powers;

    PolyDom PD(F);
    FieldPoly FP_minP = FieldPoly(minP.begin(), minP.end());
    PD.factor(factors, powers, FP_minP);

    // Factorized minP checks

    FieldPoly res;
    size_t nb_factors = factors.size();
    for(size_t i = 0; i < nb_factors; ++i)
    {
        Element_ptr E_min = FFLAS::fflas_new(F, n);
        FFLAS::fzero(F, n, E_min, 1);
        PD.div(res, FP_minP, factors[i]);

        for(size_t j = 0; j < res.size(); ++j)
            FFLAS::faxpy(F, n, res[j], K+j*ldk, 1, E_min, 1);
        if(FFLAS::fiszero(F, n, E_min, 1))
        {
            cout<<"NONMINIMALERROR"<<endl;
            return false;
        }
        FFLAS::fflas_delete(E_min);
    }


    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(K);
    return true;
}

template <class Field>
bool run_with_field (Givaro::Integer q, size_t b, size_t n, size_t iters, uint64_t seed)
{
    bool ok = true;
    int nbiter = iters;

    while (ok && nbiter)
    {
        Field* F = chooseField<Field>(q, b, seed); // F, characteristic q of b bits
        typename Field::RandIter G(*F, seed++); //random generator over F

        if(F == nullptr)
            return true; //if F is null, nothing to test, just pass

        ostringstream oss;
        F->write(oss);
        cout.fill('.');
        cout<<"Checking ";
        cout.width(40);
        cout<<oss.str();
        cout<<" ... ";

        ok = ok && check_minpoly(*F, n, G);

        if(!ok)
            cout<<"FAILED"<<endl;
        else
            cout<<"PASS"<<endl;

        delete F;
        nbiter--;
    }


    return ok;
}


int main(int argc, char** argv)
{
    /* Test parameters */
    Givaro::Integer q = -1;
    size_t b = 0;
    size_t n = 108;
    size_t iters = 1;
    bool loop = false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",	TYPE_INTEGER, &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT, &b },
        { 'n', "-n N", "Set the order of the matrix.", TYPE_INT, &n },
        { 'i', "-i, R", "set the number of repetitions.", TYPE_INT, &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "set seed for the random generator.", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc, argv, as);

    bool ok = true;

    do
    {
        ok = ok && run_with_field<Modular<double>>(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<float>>(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double>>(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float>>(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int32_t>>(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t>>(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int64_t>>(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t>>(q,b,n,iters,seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,n/6,iters,seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),n/6,iters,seed);
    } while(ok && loop);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
