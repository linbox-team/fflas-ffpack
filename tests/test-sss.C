/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet and Quentin Houssier
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

template<class Field>
bool test_diag_product (const Field & F, size_t n, size_t s, size_t t,
			typename Field::ConstElement_ptr D, size_t ldd,typename Field::Element_ptr TS, size_t ldt)
{
    typedef typename Field::Element_ptr Element_ptr ;
    Element_ptr null = fflas_new (F, n, s);
    fzero (F, n, s, null, s);
    Element_ptr C = fflas_new (F, n, t);
    
    productSSSxTS (F, n, s, null, s, null, s, null, s, null, s, null, s, null, s, D, ldd, t, F.one, TS, ldt, F.zero, C, t);

    std::cout<<"D:"<<std::endl;
    for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < s; column++)
	  {
	    std::cout<<D[column + line*s]<<" ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<"B:"<<std::endl;
    for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < t; column++)
	  {
	    std::cout<<TS[column + line*ldt]<<" ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<"C:"<<std::endl;
    for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < t; column++)
	  {
	    std::cout<<C[column + line*t]<<" ";
	  }
	std::cout<<std::endl;
      }
    FFLAS::fflas_delete(null);
    FFLAS::fflas_delete(C);
    return true;
}

/* @brief Tests SSSxTS product on a lower triangular matrix */
template<class Field>
bool test_LT_product (const Field & F, size_t n, size_t s, size_t t,
		      typename Field::ConstElement_ptr P, size_t ldp,
		      typename Field::ConstElement_ptr Q, size_t ldq,
		      typename Field::ConstElement_ptr R, size_t ldr,
		      typename Field::Element_ptr TS, size_t ldt)
{
    typedef typename Field::Element_ptr Element_ptr ;
    Element_ptr null = fflas_new (F, n, s);
    fzero (F, n, s, null, s);
    Element_ptr C = fflas_new (F, n, t);
    
    productSSSxTS (F, n, s, P, ldp, Q, ldq, R, ldr, null, s, null, s, null, s, null, s, t, F.one, TS, ldt, F.zero, C, t);

    std::cout<<"P:"<<std::endl;
    for (size_t line = 0; line < n - s; line++)
      {
	for (size_t column = 0; column < s; column++)
	  {
	    std::cout<<P[column + line*ldp]<<" ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<"Q:"<<std::endl;
    for (size_t line = 0; line < n - s; line++)
      {
	for (size_t column = 0; column < s; column++)
	  {
	    std::cout<<Q[column + line*ldq]<<" ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<"R:"<<std::endl;
    for (size_t line = 0; (line < n - 2*s) && (2*s <= n); line++)
      {
	for (size_t column = 0; column < s; column++)
	  {
	    std::cout<<R[column + line*ldq]<<" ";
	  }
	std::cout<<std::endl;
      }
    
    std::cout<<"B:"<<std::endl;
    for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < t; column++)
	  {
	    std::cout<<TS[column + line*ldt]<<" ";
	  }
	std::cout<<std::endl;
      }
    std::cout<<"C:"<<std::endl;
    for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < t; column++)
	  {
	    std::cout<<C[column + line*t]<<" ";
	  }
	std::cout<<std::endl;
      }
    
    /* Free memory and return */
    FFLAS::fflas_delete(null);
    FFLAS::fflas_delete(C);
    return true;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t s, size_t t, uint64_t seed){
  bool ok = true ;
  typedef typename Field::Element_ptr Element_ptr ;
  Field* F= chooseField<Field>(q,b,seed);
  if (F==nullptr)
    return true;
  
  Element_ptr D = fflas_new (*F, n, s);
  Element_ptr P = fflas_new (*F, n - s, s);
  Element_ptr Q = fflas_new (*F, n, s);
  Element_ptr R = fflas_new (*F, n - s, s);
  Element_ptr B = fflas_new (*F, n, t);
  
  for (size_t line = 0; line < n; line++)
      {
	for (size_t column = 0; column < s; column++)
	  {
	    D[column + line*s] = line;
	    Q[column + line * s] = 1;
	    if (line < n - s)
	      {
		P[column + line * s] = 1;
		R[column + line * s] = 2;
	      }
	  }
	for (size_t column = 0; column < t; column++)
	  {
	    B[column + line*t] = column + line;
	  }
      }
  /* Test with only Diagonal blocks */
  //  ok = ok && test_diag_product (*F, n, s, t, D, s, B, t);
  /* Test with only lower-triangular part */
  ok = ok && test_LT_product (*F, n, s, t, P, s, Q, s, R, s, B, t);

  if ( !ok )
    std::cout << "FAILED "<<std::endl;
  else
    std::cout << "PASSED "<<std::endl;
  delete F;

  /* Free memory and return */
  FFLAS::fflas_delete(D);
  FFLAS::fflas_delete(B);
  FFLAS::fflas_delete(P);
  FFLAS::fflas_delete(Q);
  FFLAS::fflas_delete(R);

  return ok;
}

int main(int argc, char** argv)
{
  cerr<<setprecision(20);
    Givaro::Integer q=101;
    bool ok=true;
    
    
    /* Very easy */
#if 0
    size_t n = 5;
    size_t s = 1;
    size_t t = 1;
#endif
    /* A bit less easy */
    size_t s = 2;
    size_t t = 2;
    size_t n = 9;

    uint64_t b = 1;

    uint64_t seed = getSeed();
    srand(seed);
    
    ok = ok && run_with_field<Givaro::Modular<int64_t> > (q, b, n, s, t, seed);
    return !ok;
}
    
#if 0

    
    fassign (F, n, n, A, lda, B, lda);

    size_t * P = fflas_new<size_t> (n);
    size_t * Q = fflas_new<size_t> (n);

    size_t r2;

    r2 =  LTBruhatGen (F, diag, n, B, lda, P, Q);
        
    size_t s = LTQSorder (n,r, P, Q);

    if (s != t){
      fail = true;
      std::cerr<<"ERROR: wrong quasi-separable order (expected "<<t<<" but got "<<s<<")"<<std::endl;
    }
    if (r2 != r){
      fail=true;
      std::cerr<<"ERROR: wrong rank (expected "<<r<<" but got "<<r2<<")"<<std::endl;
    }

    Element_ptr L = fflas_new(F,n,n);
    Element_ptr R = fflas_new(F,n,n);
    Element_ptr U = fflas_new(F,n,n);

    // TODO: later on, don't build a dense matrix for R
    getLTBruhatGen(F, n, r, P, Q, R, n);
    getLTBruhatGen(F, FflasLower, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, n, r, P, Q, B, lda, L,n);
    getLTBruhatGen(F, FflasUpper, diag, n, r, P, Q, B, lda, U, n);

     // WriteMatrix(std::cerr<<"L = "<<std::endl,F,n,n,L,n);
     // WriteMatrix(std::cerr<<"U = "<<std::endl,F,n,n,U,n);

        //test of compression into block bi diagonal
    Element_ptr U2= fflas_new(F, n,n);
    Element_ptr Xu = fflas_new(F, 2*s, n);
    size_t * Ku = fflas_new<size_t> (r+1);
    size_t * Mu = fflas_new<size_t> (n);
    size_t * Tu = fflas_new<size_t>(r);
    size_t NbBlocksU = CompressToBlockBiDiagonal(F, FflasUpper, n, s, r, P, Q, U,n ,Xu,n,Ku,Mu,Tu);
    ExpandBlockBiDiagonalToBruhat(F,FflasUpper,n,s,r,U2,n,Xu,n,NbBlocksU,Ku,Mu,Tu);
    if(!fequal(F,n,n,U,n,U2,n)){
      fail= true;
      std::cerr<<"ERROR: Compression of U lost information"<<std::endl;
    }
    Element_ptr L2= fflas_new(F, n,n);
    Element_ptr Xl = fflas_new(F, n, 2*s);
    size_t * Kl = fflas_new<size_t> (r+1);
    size_t * Ml = fflas_new<size_t> (n);
    size_t * Tl = fflas_new<size_t>(r);
    size_t NbBlocksL = CompressToBlockBiDiagonal(F, FflasLower, n, s, r, P, Q, L,n ,Xl,2*s,Kl,Ml,Tl);
    ExpandBlockBiDiagonalToBruhat(F,FflasLower,n,s,r,L2,n,Xl,2*s,NbBlocksL,Kl,Ml,Tl);
    if(!fequal(F,n,n,L,n,L2,n)){
      fail= true;
      std::cerr<<"ERROR: Compression of L lost information"<<std::endl;
    }
    // WriteMatrix(std::cerr<<"A="<<std::endl,F,n,n, A,lda)<<std::endl;
    // WriteMatrix(std::cerr<<"TS="<<std::endl,F,n,l, TS,l)<<std::endl;

    Element_ptr CBruhat = fflas_new(F, n, l);
    productBruhatxTS(F, n, s, r, l, P, Q, Xu, n, NbBlocksU, Ku, Tu, Mu,Xl, 2*s, NbBlocksL, Kl, Tl, Ml,TS, l, F.zero, CBruhat, l);
    Element_ptr Cfgemm = fflas_new(F, n, l);
    fgemm(F, FflasNoTrans, FflasNoTrans, n,l,n,F.one,A, lda, TS, l, F.zero, Cfgemm, l); 

    if(!fequal(F,n,l,CBruhat,l,Cfgemm,l)){
      fail= true;
      std::cerr<<"ERROR: fgmemm != productBruhatxTS"<<std::endl;
      WriteMatrix(std::cerr<<"CBruhat = "<<std::endl,F,n,l,CBruhat,l)<<std::endl;
      WriteMatrix(std::cerr<<"Cfgemm = "<<std::endl,F,n,l,Cfgemm,l)<<std::endl;
    }

        //fflas_delete ( U2, Xu,Ku,Mu,Tu,L2,Kl,Xl,Ml,Tl);
    // B <- L R^T
    fgemm(F, FflasNoTrans, FflasTrans, n,n,n, F.one, L, n, R, n, F.zero, B, lda);
    // L <- B U
    fgemm(F, FflasNoTrans, FflasNoTrans, n,n,n, F.one, B, lda, U, n, F.zero, L, n);
    // Extract the left triangular part of L
    for (size_t i=0; i<n; ++i)
      fzero(F, i+1, L + i*n + n-i-1, 1);
    // Check L == A
    if (!fequal (F, n, n, L, n, A, lda)){
      fail = true;
      std::cerr<<"ERROR: A != Left(L R^T U)"<<std::endl;
    }

    fflas_delete (B, P, Q);
    return fail;
}

template<class Field, FFLAS_DIAG diag, class RandIter>
bool launch_test (const Field & F, size_t n, size_t r, size_t t, size_t l,RandIter& G)
{
    //typedef typename Field::Element Element ;
    typedef typename Field::Element_ptr Element_ptr ;
    bool fail = false ;
    { /*  user provided params, larger lda */
        size_t lda = n+10 ;
        Element_ptr A = fflas_new (F, n, lda);
	Element_ptr TS = fflas_new(F, n, l);
            // TODO implement this randomGenerator
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A,lda,G);
	RandomMatrix(F, n, l, TS,l,G);
        fail = fail || test_BruhatGenerator <Field,diag> (F, n, r, t, A, lda,TS, l,G);

        if (fail) std::cout << "failed at user params" << std::endl;
        fflas_delete( A );
    }
    { /*  user provided n, larger lda, large small t */
    }
    { /*  user provided n, larger lda, large large t */
    }
    return !fail;
}
template <class Field, class RandGen>
bool testLTQSRPM (const Field & F,size_t n, size_t r, size_t t, RandGen& G){

    size_t * rows = FFLAS::fflas_new<size_t>(r);
    size_t * cols = FFLAS::fflas_new<size_t>(r);
    RandomLTQSRankProfileMatrix (n, r,  t, rows, cols);

    typename Field::Element_ptr A = fflas_new(F,n,n);
    getLTBruhatGen(F, n, r, rows, cols, A, n);

        //WriteMatrix (std::cerr<<"A = "<<std::endl,F,n,n,A,n);
    
    fflas_delete(A);
    size_t s = LTQSorder (n, r, rows, cols);
    if (s==t){
//        std::cerr<<"PASS"<<std::endl;
        return true;
    } else {
        std::cerr<<"Failed testLTQSRPM: QS order expected: "<<t<<", but got "<<s<<std::endl;
        return false;
    }
}




    size_t b=0;
    size_t n=93;
    size_t r=30;
    size_t t=8;
    size_t m=6;
    size_t iters=3;
    bool loop=false;


    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'n', "-n N", "Set the matrix order.", TYPE_INT , &n },
        { 'r', "-r R", "Set the rank.", TYPE_INT , &r },
        { 'm', "-m M", "Set the col dim of the TS matrix.", TYPE_INT , &m },
        { 't', "-t T", "Set the order of quasi-separability.", TYPE_INT , &t },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    if (r > n) r = n/2;
    if (t > r) t = r/2;


    

    do{
        ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,r,t,m,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,n/4,r/4,t/4,m/4,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),n/4,r/4,t/4,m/4,iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;
    
    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

#endif
