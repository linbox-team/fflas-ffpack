/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Hippolyte Signargout (hippolyte.signargout@ens-lyon.fr)
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


//-------------------------------------------------------------------------
//      Test suite for the Quasi-Separable matrices in SSS format
//-------------------------------------------------------------------------

/* Structure taken from test-quasisep.C */

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

#include <random>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

/** \brief test equality between matrix reconstructed by SSSToDense and applied to identity with productSSSxTS */
template<class Field>
bool test_reconstruction_compatibility (const Field & F, size_t n, size_t s, 
                                        typename Field::ConstElement_ptr P, size_t ldp,
                                        typename Field::ConstElement_ptr Q, size_t ldq,
                                        typename Field::ConstElement_ptr R, size_t ldr,
                                        typename Field::ConstElement_ptr U, size_t ldu,
                                        typename Field::ConstElement_ptr V, size_t ldv,
                                        typename Field::ConstElement_ptr W, size_t ldw,
                                        typename Field::ConstElement_ptr D, size_t ldd)
{
    typename Field::Element_ptr rec = fflas_new (F, n, n);
    typename Field::Element_ptr app = fflas_new (F, n, n);
    typename Field::Element_ptr Id = fflas_new (F, n, n);
    fidentity (F, n, n, Id, n);
  
    SSSToDense (F, n, s, P, ldp, Q, ldq, R, ldr, U, ldu, V, ldv, W, ldw,
                D, ldd, rec, n);
    productSSSxTS (F, n, n, s, F.one, P, ldp, Q, ldq, R, ldr, U, ldu, V, ldv, W, ldw, D, ldd,
                   Id, n, F.zero, app, n);

    bool ok = fequal (F, n, n, rec, n, app, n);

    if ( !ok )
    {
        std::cout << "ERROR: different results for reconstruction and application to identity "<<std::endl;
        WriteMatrix(std::cout<<"rec = "<<std::endl, F, n, n, rec, n);
        WriteMatrix(std::cout << "app =  "<<std::endl, F, n, n, app, n);
    }

    FFLAS::fflas_delete(rec);
    FFLAS::fflas_delete(app);
    FFLAS::fflas_delete(Id);

    return ok;
}

/** \brief test equality between applying productSSSxTS or dense matrix */
template<class Field>
bool test_application_compatibility (const Field & F, size_t n, size_t t, size_t s,
                                     const typename Field::Element alpha,
                                     typename Field::ConstElement_ptr P, size_t ldp,
                                     typename Field::ConstElement_ptr Q, size_t ldq,
                                     typename Field::ConstElement_ptr R, size_t ldr,
                                     typename Field::ConstElement_ptr U, size_t ldu,
                                     typename Field::ConstElement_ptr V, size_t ldv,
                                     typename Field::ConstElement_ptr W, size_t ldw,
                                     typename Field::ConstElement_ptr D, size_t ldd,
                                     typename Field::ConstElement_ptr B, size_t ldb,
                                     const typename Field::Element beta,
                                     typename Field::ConstElement_ptr C, size_t ldc)
{
    typename  Field::Element_ptr A = fflas_new (F, n, n);
    SSSToDense (F, n, s, P, ldp, Q, ldq, R, ldr, U, ldu, V, ldv, W, ldw,
                D, ldd, A, n);
  
    typename  Field::Element_ptr dense = fflas_new (F, n, t);
    typename  Field::Element_ptr qs = fflas_new (F, n, t);
    fassign (F, n, t, C, ldc, dense, t);
    fassign (F, n, t, C, ldc, qs, t);

    fgemm (F, FflasNoTrans, FflasNoTrans, n, t, n, alpha, A, n, B, ldb, beta, dense, t);
    productSSSxTS (F, n, t, s, alpha, P, ldp, Q, ldq, R, ldr, U, ldu, V, ldv, W, ldw, D, ldd,
                   B, ldb, beta, qs, t);

    bool ok = fequal (F, n, t, dense, t, qs, t);
 if ( !ok )
        std::cout << "ERROR: different results for dense and qs application "<<std::endl;
   
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(dense);
    FFLAS::fflas_delete(qs);

    return ok;
}
/** \brief test equality between block diagonal and result of densifying and 
    compressing it */
template<class Field>
//bool test_diagonal_compression (const Field & F, size_t n, size_t s,
bool test_compression (const Field & F, size_t n, size_t s,
                                     typename Field::ConstElement_ptr P, size_t ldp,
                                     typename Field::ConstElement_ptr Q, size_t ldq,
                                     typename Field::ConstElement_ptr R, size_t ldr,
                             typename Field::ConstElement_ptr U, size_t ldu,
                             typename Field::ConstElement_ptr V, size_t ldv,
                             typename Field::ConstElement_ptr W, size_t ldw,
                             typename Field::ConstElement_ptr D, size_t ldd)
{
    size_t rs = n%s;           // Size of the partial block
    size_t ls = (rs)? rs: s;   // Size of the last block
    typename  Field::Element_ptr A1 = fflas_new (F, n, n);
    typename  Field::Element_ptr A2 = fflas_new (F, n, n);
    //    typename  Field::Element_ptr Z = fflas_new (F, n + s, s); //To be sure
    typename  Field::Element_ptr Pcheck = fflas_new (F, n - s, s);
    typename  Field::Element_ptr Qcheck = fflas_new (F, n - ls, s);
    typename  Field::Element_ptr Rcheck = fflas_new (F, ((n > s + ls)? (n - s - ls): 0), s);
    typename  Field::Element_ptr Dcheck = fflas_new (F, n, s);
    typename  Field::Element_ptr Ucheck = fflas_new (F, n - ls, s);
    typename  Field::Element_ptr Vcheck = fflas_new (F, n - ls, s);
    typename  Field::Element_ptr Wcheck = fflas_new (F, ((n > s + ls)? (n - s - ls): 0), s);
    //    fzero (F, n + s, s, Z, s);

    SSSToDense (F, n, s, P, ldp, Q, ldq, R, ldq, U, ldu, V, ldv, W, ldw,
                D, ldd, A1, n);
    DenseToSSS (F, n, s, Pcheck, s, Qcheck, s, Rcheck, s, Ucheck, s, Vcheck, s, Wcheck, s,
                Dcheck, s, A1, n);
    SSSToDense (F, n, s, Pcheck, s, Qcheck, s, Rcheck, s, Ucheck, s, Vcheck, s, Wcheck, s,
                Dcheck, s, A2, n);

    bool ok = fequal (F, n, n, A1, n, A2, n);
 if ( !ok )
        {
        std::cout << "ERROR: different results for dense from generator and compression"
                  <<std::endl;
        WriteMatrix(std::cout<<"A1 = "<<std::endl, F, n, n, A1, n);
        WriteMatrix(std::cout << "A2 =  "<<std::endl, F, n, n, A2, n);
        WriteMatrix(std::cout << "Ucheck = "<<std::endl, F, n - ls, s, Ucheck, s);
        WriteMatrix(std::cout << "Vcheck =  "<<std::endl, F, n - ls, s, Vcheck, s);
        WriteMatrix(std::cout << "Wcheck = " <<std::endl, F, ((n > s + ls)? (n - s - ls): 0), s, Wcheck, s);
        WriteMatrix(std::cout << "Dcheck = " <<std::endl, F, n, s, Dcheck, s);
        }
   
 FFLAS::fflas_delete(A1, A2, Dcheck, Ucheck, Vcheck, Wcheck, Pcheck, Qcheck, Rcheck);
    return ok;
}

template<class Field>
bool launch_instance_check (const Field& F, size_t n, size_t s, size_t t, typename Field::RandIter& G)
{
        /* Generate generators */
    size_t rs = n%s;           // Size of the partial block
    size_t ls = (rs)? rs: s;   // Size of the last block
    //std::cout << "n = " << n << std::endl << "s = " << s << std::endl << "ls = " << ls << std::endl;
    typedef typename Field::Element_ptr Element_ptr;
    Element_ptr D = fflas_new (F, n, s);
    Element_ptr P = fflas_new (F, n - s, s);
    Element_ptr Q = fflas_new (F, n - ls, s);
    Element_ptr R = fflas_new (F, ((n > s + ls)? (n - s - ls): 0), s);
    Element_ptr U = fflas_new (F, n - ls, s);
    Element_ptr V = fflas_new (F, n - ls, s);
    Element_ptr W = fflas_new (F, ((n > s + ls)? (n - s - ls): 0), s);
    Element_ptr C = fflas_new (F, n, t);
    Element_ptr B = fflas_new (F, n, t);

    frand (F, G, n, s, D, s);
    frand (F, G, n - s, s, P, s);
    frand (F, G, n - ls, s, Q, s);
    frand (F, G, ((n > s + ls)? (n - s - ls): 0), s, R, s);
    frand (F, G, n - ls, s, U, s);
    frand (F, G, n - ls, s, V, s);
    frand (F, G, ((n > s + ls)? (n - s - ls): 0), s, W, s);
    frand (F, G, n, t, C, t);
    frand (F, G, n, t, B, t);

    typename Field::Element alpha, beta;
    G.random(alpha);
    G.random(beta);

    bool ok = true;
        /* Call to functions being implemented */
    ok = ok && test_reconstruction_compatibility(F, n, s, P, s, Q, s, R, s, U, s, V, s, W, s, D, s);
    ok = ok && test_application_compatibility(F, n, t, s, alpha, P, s, Q, s, R, s, U, s, V, s, W, s, D, s, B, t, beta, C, t);
    ok = ok && test_application_compatibility(F, n, t, s, alpha, P, s, Q, s, R, s, U, s, V, s, W, s, D, s, B, t, F.zero, C, t);
    ok = ok && test_compression (F, n, s, P, s, Q, s, R, s, U, s, V, s, W, s, D, s);
    
    if ( !ok )
    {
        std::cout << "FAILED "<<std::endl;
            /* Print generators for debugging */
        WriteMatrix (std::cout << "P =  "<<std::endl,F, n - ls, s, P, s);
        WriteMatrix (std::cout << "Q =  "<<std::endl,F, n - ls, s, Q, s);
        WriteMatrix (std::cout << "R =  "<<std::endl,F, n - s - ls, s, R, s);
        WriteMatrix (std::cout << "U =  "<<std::endl,F, n - ls, s, U, s);
        WriteMatrix (std::cout << "V =  "<<std::endl,F, n - ls, s, V, s);
        WriteMatrix (std::cout << "W =  "<<std::endl,F, n - s - ls, s, W, s);
        WriteMatrix (std::cout << "D =  "<<std::endl,F, n, s, D, s);
    }
        /* Free memory and return */
    FFLAS::fflas_delete(D);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(P);
    FFLAS::fflas_delete(Q);
    FFLAS::fflas_delete(R);
    FFLAS::fflas_delete(C);
    FFLAS::fflas_delete(U);
    FFLAS::fflas_delete(V);
    FFLAS::fflas_delete(W);

    return ok;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t s, size_t t, size_t iters, uint64_t seed)
{
    bool ok = true ;
    while (ok && iters)
    {
            //std::cout << "Iterations left: " << iters << std::endl;
            /* New field 
             * chooseField returns a pointer, F needs to be passed by its value */
        Field* F= chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
            /* Initiate random number generator */
        typename Field::RandIter G(*F,seed++);
        srandom(seed);
        
        std::ostringstream oss;
        F->write(oss);

        std::cout.fill('.');
        std::cout<<"Checking ";
        std::cout.width(65);
        std::cout<<oss.str();
        std::cout<<" ... ";

        ok = ok && launch_instance_check (*F, n, s, t, G);
        ok = ok && launch_instance_check (*F, n, (random() % s)+1, t, G);

        if (ok)
            std::cout << "PASSED "<<std::endl;

        delete F;
        iters--;
    }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(20); // In order to print integers as integers even on float types, could be done once for all fflas
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=319;
    size_t s=23;
    size_t t=42;
    int iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT , &b },
        { 'n', "-n N", "Set the matrix row and column dimension.", TYPE_INT , &n },
        { 's', "-s S", "Set the order of quasi-separability.", TYPE_INT , &s },
        { 't', "-t T", "Set the col dim of the Tall and Skinny matrix.", TYPE_INT , &t },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT, &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
   
    bool ok=true;
    do{
        ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,n,s,t,iters,
                                                                     seed); // Valgrind does not like this one 
        ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,s,t,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,9, ceil(n/4.), ceil(s / 4.), ceil(t / 4.), iters,
                                                                     seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:224), ceil(n/4.), ceil(s / 4.), ceil(t / 4.),
                                                                     iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s