/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (c) FFLAS-FFPACK
* Written by Bastien Vialla <bastien.vialla@lirmm.fr>
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
*/

#define __DLP_CHALLENGE

#include <iostream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include "gmpxx.h"
#include <givaro/zring.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/givinteger.h>
#include <recint/recint.h>
#include <givaro/givintprime.h>

using namespace RecInt;

#include "fflas-ffpack/fflas/fflas_sparse.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/fflas/fflas_sparse/read_sparse.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/flimits.h"

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif

using namespace std;
using namespace FFLAS;
using namespace Givaro;

/*******************************************************************************************************************
 *
 *      Utility functions: sms reader and random field
 *  
 *******************************************************************************************************************/

template<class T>
size_t bitSize(T n){
  return sizeof(T)*4-__builtin_clz(n);
}

template<typename Field>
Givaro::Integer maxFieldElt() {return (Givaro::Integer)Field::maxCardinality();} 
template<>
Givaro::Integer maxFieldElt<Givaro::ZRing<Givaro::Integer>>() {return (Givaro::Integer)-1;}

/*** Field chooser for test according to characteristic q and bitsize b ***/
/* if q=-1 -> field is chosen randomly with a charateristic of b bits
   if b=0 -> bitsize is chosen randomly according to maxFieldElt
*/
template<typename Field>
Field* chooseField(Givaro::Integer q, uint64_t b){
  Givaro::Integer maxV= maxFieldElt<Field>();
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 mt_rand(seed);
  if (maxV>0 && (q> maxV || b> maxV.bitsize()))
    return nullptr;
  if (b<=1){
    //srand((double)std::chrono::high_resolution_clock::now());
    auto bitrand = std::bind(std::uniform_int_distribution<uint64_t>(2,maxV.bitsize()-1),
                 mt_rand);
    b = bitrand();
  }
  Givaro::IntPrimeDom IPD;
  Givaro::Integer tmp,p;
  if (q==-1){
    // Choose characteristic as a random prime of b bits
    do{
      Givaro::Integer _p;
      Givaro::Integer::seeding(Givaro::Integer(mt_rand()));
      Givaro::Integer::random_exact_2exp(_p,b);
      IPD.prevprime( tmp, _p+1 );
      p =  tmp;
    }while( (p < 2) );
  }
  else p=q;

  return new Field(p);
}

/*************************************************************************************************************/

int main(int argc, char **argv) {
    using Field        = Modular<Integer>;
    using FieldMat     = ZRing<double>;
    using FieldComp    = FFPACK::RNSIntegerMod<FFPACK::rns_double>;
    using FieldElement = RecInt::rmint<7>;
    using FieldRec     = ZRing<FieldElement>;
    using SparseMatrix = FFLAS::Sparse<FieldRec, FFLAS::SparseMatrix_t::HYB_ZO>;

    Integer q = -1;
    int b = 128;
    int blockSize = 1;
    std::string matrixFile = "";
    int nIter = 100;
    
    static Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",   TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",   TYPE_INT , &b },
        { 'k', "-k K", "Set the size of the block (1 by default).",       TYPE_INT, &blockSize },
        { 'n', "-n N", "Set the size of the block (1 by default).",       TYPE_INT, &nIter },
        { 'f', "-f FILE", "Set matrix file.",                             TYPE_STR, &matrixFile },
         END_OF_ARGUMENTS };

    FFLAS::parseArguments(argc, argv, as);

    // Construct Givaro::Integer field
    Field *F= chooseField<Field>(q,b);
    if (F==nullptr) exit(0);
    Integer p;
    F->cardinality(p);
    cout << "Prime p: " << p << endl;

    RecInt::ruint<7> pRec;
    // RecInt::mpz_to_ruint(pRec, FieldElement(p));
    FieldElement::init_module(ruint<7>(p));
    FieldRec Frec;
    // Pointers for the matrix
    index_t *row = nullptr, *col = nullptr;
    typename FieldRec::Element_ptr dat;
    index_t rowdim, coldim;
    uint64_t nnz;
  
    // Read the matrix
    readSmsFormat(matrixFile, Frec, row, col, dat, rowdim, coldim, nnz);
    vector<index_t> rowCoo(nnz, 0);
    for(size_t i = 0 ; i < rowdim ; ++i){
        for(size_t j = row[i] ; j < row[i+1] ; ++j){
            rowCoo[j] = i;
        }
    }

    // Build the matrix
    SparseMatrix A;
    FFLAS::sparse_init(Frec, A, rowCoo.data(), col, dat, rowdim, coldim, nnz);

    FFLAS::fflas_delete(row);
    FFLAS::fflas_delete(col);
    FFLAS::fflas_delete(dat);
    rowCoo.resize(0);

    vector<FieldElement> x(coldim*blockSize, 1), y(rowdim*blockSize, 0);

    pfspmm(Frec, A, blockSize, x.data(), blockSize, 0, y.data(), blockSize);

    return 0;
}

