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

#include <iostream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include "gmpxx.h"
#include "givaro/unparametric.h"
#include "givaro/modular.h"
#include "givaro/modular-balanced.h"
#include "givaro/givinteger.h"
#include "recint/recint.h"
#include <givaro/givintprime.h>

#include "fflas-ffpack/fflas/fflas_sparse.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/field/rns-integer-mod.h"
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

/*
 * Reader for sms format.
 * Return 3 arrays of size nnz:
 *  - val: values of non zeros elements
 *  - col: column index of non zeros element
 *  - row: row index of non zeros element
 */
template <class Field>
void readSmsFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
                   typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz, bool trans = false, bool unparam = false) {
    using namespace details_spmv;
    std::ifstream file(path, std::ios::in);
    std::vector<std::string> tokens;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::copy(std::istream_iterator<std::string>(is), std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string>>(tokens));
    if(trans){
        coldim = static_cast<index_t>(std::stoull(tokens[0]));
        rowdim = static_cast<index_t>(std::stoull(tokens[1]));
    }else{
        rowdim = static_cast<index_t>(std::stoull(tokens[0]));
        coldim = static_cast<index_t>(std::stoull(tokens[1]));
    }
    
    std::vector<Coo<Field>> data;
    nnz = 0;
    while (std::getline(file, line)) {
        tokens.resize(0);
        std::istringstream iss(line);

        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string>>(tokens));

        if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
            typename Field::Element v;
            if(unparam){
                v = std::stol(tokens[2]);
            }else{
                f.init(v, std::stol(tokens[2]));
            }
            index_t r = (index_t)(std::stoull(tokens[0])) - 1;
            index_t c = (index_t)(std::stoull(tokens[1])) - 1;
            if(trans){
                data.emplace_back(v, c, r);
            }else{
                data.emplace_back(v, r, c);
            }
        }
    }
    std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });

    row = FFLAS::fflas_new<index_t>(data.size());
    col = FFLAS::fflas_new<index_t>(data.size());
    val = FFLAS::fflas_new(f, data.size(), 1);
    nnz = data.size();

    for (size_t i = 0, end = data.size(); i < end; ++i) {
        val[i] = data[i].val;
        col[i] = data[i].col;
        row[i] = data[i].row;
    }
}

template<class T>
size_t bitSize(T n){
  return sizeof(T)*4-__builtin_clz(n);
}

template<typename Field>
Givaro::Integer maxFieldElt() {return (Givaro::Integer)Field::getMaxModulus();} 
template<>
Givaro::Integer maxFieldElt<Givaro::UnparametricRing<Givaro::Integer>>() {return (Givaro::Integer)-1;}

/*** Field chooser for test according to characteristic q and bitsize b ***/
/* if q=-1 -> field is chosen randomly with a charateristic of b bits
   if b=0 -> bitsize is chosen randomly according to maxFieldElt
*/
template<typename Field>
Field* chooseField(Givaro::Integer q, unsigned long b){
  Givaro::Integer maxV= maxFieldElt<Field>();
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 mt_rand(seed);
  if (maxV>0 && (q> maxV || b> maxV.bitsize()))
    return nullptr;
  if (b<=1){
    //srand((double)std::chrono::high_resolution_clock::now());
    auto bitrand = std::bind(std::uniform_int_distribution<unsigned long>(2,maxV.bitsize()-1),
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
    using FieldMat     = UnparametricRing<double>;
    using FieldComp    = FFPACK::RNSIntegerMod<FFPACK::rns_double>;
    using SparseMatrix = FFLAS::Sparse<FieldMat, FFLAS::SparseMatrix_t::CSR_HYB>;

    size_t iter = 10;
    Integer q = -1;
    int b = 128;
    int blockSize = 1;
    std::string matrixFile = "";
    
    static Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",   TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",   TYPE_INT , &b },
        { 'k', "-k K", "Set the size of the block (1 by default).",       TYPE_INT, &blockSize },
        { 'f', "-f FILE", "Set matrix file.",                             TYPE_STR, &matrixFile },
         END_OF_ARGUMENTS };

    FFLAS::parseArguments(argc, argv, as);

    // Construct Givaro::Integer field
    Field *F= chooseField<Field>(q,b);
    if (F==nullptr) exit(0);
    Integer p;
    F->cardinality(p);
  
    // Pointers for the matrix
    index_t *row = nullptr, *col = nullptr;
    typename FieldMat::Element_ptr dat;
    index_t rowdim, coldim;
    uint64_t nnz;
    // Field associate to the matrix
    FieldMat Fword;
  
    // Read the matrix
    readSmsFormat(matrixFile, Fword, row, col, dat, rowdim, coldim, nnz, true, false);
    
    // Build the matrix
    SparseMatrix A;
    FFLAS::sparse_init(Fword, A, row, col, dat, rowdim, coldim, nnz);

    FFLAS::fflas_delete(row);
    FFLAS::fflas_delete(col);
    FFLAS::fflas_delete(dat);

    vector<double> x(coldim, 1), y(rowdim, 0);

    // Compute the bigger row
    FFLAS::fspmv(Fword, A, x.data(), 0, y.data());
    double maxSum = *(std::max_element(y.begin(), y.end()));
    cout << "maxSum: " << maxSum << endl;

    // Compute the bitsize of the RNS primes
    size_t primeBitsize = 53 - Integer(maxSum).bitsize();
    cout << "primeBitsize: " << primeBitsize << endl;
    // construct RNS
    FFPACK::rns_double RNS(Integer(maxSum)*p, primeBitsize,true);
    size_t rnsSize = RNS._size;
    
    // construct RNS field
    FieldComp Frns(p,RNS);

    std::vector<Integer> X(coldim*blockSize), Y(rowdim*blockSize, 0);

    // Fill X with random values
    for(auto &x: X){
        Givaro::Integer::random_exact_2exp(x,b);
        F->init(x, x);
    }

    double* Xrns = fflas_new<double>(coldim*blockSize*rnsSize, Alignment::CACHE_LINE);
    double* Yrns = fflas_new<double>(rowdim*blockSize*rnsSize, Alignment::CACHE_LINE);

    // Transform X in RNS
    RNS.init_dlp(coldim*blockSize, Xrns, X.data(), 1);

    // perform Yrns = A.Xrns + beta.Yrns over ZZ
    fspmm(Fword, A, blockSize*rnsSize, Xrns, blockSize*rnsSize, 0, Yrns, blockSize*rnsSize);

    // reduce Yrns wrt the RNS basis
    RNS.reduce(rowdim*blockSize, Yrns, 1, true);

    // reduce Yrns mod p
    Frns.reduce_modp_rnsmajor_scal_quad(rowdim*blockSize,  FFPACK::rns_double_elt_ptr(Yrns, 1));

    // Reconstruct Y from Yrns
    RNS.convert_dlp(rowdim*blockSize, Y.data(), 1, Yrns);

    FFLAS::fflas_delete(Xrns);
    FFLAS::fflas_delete(Yrns);

    return 0;
}

