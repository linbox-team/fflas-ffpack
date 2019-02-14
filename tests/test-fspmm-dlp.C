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

using Data = std::vector<details_spmv::Coo<ZRing<double>>>;
using Coo = typename Data::value_type;

/*******************************************************************************************************************
 *
 *      Utility functions: sms reader and random field
 *
 *******************************************************************************************************************/

void readMat(string path, index_t *& row, index_t *& col, double *&val, index_t &rowdim, index_t &coldim, uint64_t & nnz){
    std::ifstream file(path, std::ios::out);
    std::string line, nnz_c;
    std::getline(file, line);
    std::istringstream(line) >> rowdim >> coldim >> nnz_c;
    Data mat;
    int64_t r, c, v;
    while(std::getline(file, line)){
        std::istringstream(line) >> r >> c >> v;
        if(r!=0)
            mat.emplace_back(v, r-1,c-1);
    }
    std::sort(mat.begin(), mat.end(),
              [](Coo &a, Coo &b){
              return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
              ;});
    mat.shrink_to_fit();
    nnz = mat.size();
    val = fflas_new<double>(nnz, Alignment::CACHE_LINE);
    col = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
    row = fflas_new<index_t>(nnz, Alignment::CACHE_LINE);
    for(size_t i = 0 ; i < nnz ; ++i){
        val[i] = mat[i].val;
        col[i] = mat[i].col;
        row[i] = mat[i].row;
    }
}

/*************************************************************************************************************/

int main(int argc, char **argv) {
    using Field        = Modular<Integer>;
    using FieldMat     = ZRing<double>;
    using FieldComp    = FFPACK::RNSIntegerMod<FFPACK::rns_double_extended>;
    using SparseMatrix = Sparse<FieldMat, SparseMatrix_t::CSR>;
    uint64_t seed = getSeed();

    Integer q = -1;
    int b = 128;
    int blockSize = 1;
    std::string matrixFile = "";
    int nIter = 100;

    static Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",   TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",   TYPE_INT , &b },
        { 'k', "-k K", "Set the size of the block (1 by default).",       TYPE_INT, &blockSize },
        { 'n', "-n N", "Number of iterations (1 by default).",       TYPE_INT, &nIter },
        { 'f', "-f FILE", "Set matrix file.",                             TYPE_STR, &matrixFile },
        { 's', "-s seed", "Set seed for the random generator", TYPE_INT, &seed },
        END_OF_ARGUMENTS };

    parseArguments(argc, argv, as);

    // Construct Givaro::Integer field
    Field *F= chooseField<Field>(q,b,seed);
    if (F==nullptr) exit(0);
    Integer p;
    F->cardinality(p);
    cout << "Prime p: " << p << endl;

    // Pointers for the matrix
    index_t *row = nullptr, *col = nullptr;
    typename FieldMat::Element_ptr dat;
    index_t rowdim, coldim;
    uint64_t nnz;
    // Field associate to the matrix
    FieldMat Fword;

    // Read the matrix
    readMat(matrixFile, row, col, dat, rowdim, coldim, nnz);

    vector<int64_t> rows(rowdim, 0);
    for(size_t i = 0 ; i < nnz ; ++i)
        rows[row[i]]++;
    for(size_t i = 0 ; i < 20 ; ++i)
        cout << "#rows with "<<i<<" nnz: "  << std::count(rows.begin(), rows.end(), i) << endl;

    // Build the matrix
    SparseMatrix A;
    sparse_init(Fword, A, row, col, dat, rowdim, coldim, nnz);

    fflas_delete(row);
    fflas_delete(col);
    fflas_delete(dat);

    vector<double> x(coldim, 1), y(rowdim, 0);

    cout.precision(20);

    // Compute the bigger row
    fspmv(Fword, A, x.data(), 0, y.data());
    for(auto &x: y){
        if(x < 0){
            x = -x;
        }
    }
    double maxSum = *(std::max_element(y.begin(), y.end()));
    cout << "maxSum: " << maxSum << endl;

    // Compute the bitsize of the RNS primes
    size_t primeBitsize = 53 - Integer(maxSum).bitsize()-1;
    cout << "primeBitsize: " << primeBitsize << endl;
    // construct RNS
    // primeBitsize = 23;
    FFPACK::rns_double_extended RNS(Integer(maxSum)*p, primeBitsize, true, 0);
    size_t rnsSize = RNS._size;
    cout << "M: " << RNS._M << endl;
    cout << "RNS basis size: " << rnsSize << endl;
    cout << "Rns basis: ";
    for(auto&x:RNS._basis){
        cout << x << " ";
    }
    cout << endl;
    cout << "RNS Mi: " << endl;
    for(auto &x : RNS._Mi){
        cout << x << " ";
    }
    cout << endl;
    cout << "RNS MMi: " << endl;
    for(auto &x : RNS._MMi){
        cout << x << " ";
    }
    cout << endl;
    // construct RNS field
    FieldComp Frns(p,RNS);

    std::vector<Integer> X(coldim*blockSize), Y(rowdim*blockSize, 0);

    // Fill X with random values
    for(auto &x: X){
        Givaro::Integer::random_exact_2exp(x,b);
        F->init(x, x);
    }

    size_t ld = 0;
    Integer maxRep = Integer(maxSum)*rnsSize*p;
    while(maxRep.bitsize() < RNS._M.bitsize()){
        maxRep *= Integer(maxSum);
        ld++;
    }
    ld -= 1;
    cout << "Spmm by modp: " << ld << endl;

    double* Xrns = fflas_new<double>(coldim*blockSize*rnsSize, Alignment::CACHE_LINE);
    double* Yrns = fflas_new<double>(rowdim*blockSize*rnsSize, Alignment::CACHE_LINE);

    // Transform X in RNS
    RNS.init(coldim*blockSize, Xrns, X.data(), 1);

    cout << endl;
    TTimer Tspmm;
    TTimer Tmodp;
    TTimer Ttotal;
    double spmmTime = 0, modpTime = 0;
    bool bb = true;
    Ttotal.start();
    for(size_t kk = 1 ; kk <= nIter ; ++kk){
        // perform Yrns = A.Xrns + beta.Yrns over ZZ
        Tspmm.start();
        if(bb){
            pfspmm(Fword, A, blockSize*rnsSize, Xrns, blockSize*rnsSize, 0, Yrns, blockSize*rnsSize);
            RNS.reduce(rowdim*blockSize, Yrns, 1, true);
            // reduce Yrns wrt the RNS basis
            Tspmm.stop();
            spmmTime += Tspmm.usertime();

            cout << "after spmm:" << endl;
            for(size_t i = 0, end = (Y.size()>20)?20:Y.size() ; i < end ; ++i){
                cout << Yrns[i] << " ";
            }
            cout << endl;
            bb = !bb;
            // if(kk%ld == 0){
            Tmodp.start();
            Frns.reduce_modp_rnsmajor_scal_quad(rowdim*blockSize,  FFPACK::rns_double_elt_ptr(Yrns, 1));
            Tmodp.stop();
            modpTime += Tmodp.usertime();
            cout << "after modp:" << endl;
            for(size_t i = 0, end = (Y.size()>20)?20:Y.size() ; i < end ; ++i){
                cout << Yrns[i] << " ";
            }
            cout << endl;
            // }
        }else{
            fspmm(Fword, A, blockSize*rnsSize, Yrns, blockSize*rnsSize, 0, Xrns, blockSize*rnsSize);
            RNS.reduce(rowdim*blockSize, Xrns, 1, true);
            // reduce Yrns wrt the RNS basis
            Tspmm.stop();
            spmmTime += Tspmm.usertime();
            bb = !bb;
            for(size_t i = 0, end = (Y.size()>20)?20:Y.size() ; i < end ; ++i){
                cout << Xrns[i] << " ";
            }
            cout << endl;
            // if(kk%ld == 0){
            Tmodp.start();
            Frns.reduce_modp_rnsmajor_scal_quad(rowdim*blockSize,  FFPACK::rns_double_elt_ptr(Xrns, 1));
            Tmodp.stop();
            modpTime += Tmodp.usertime();
            // }
            cout << "after modp:" << endl;
            for(size_t i = 0, end = (Y.size()>20)?20:Y.size() ; i < end ; ++i){
                cout << Xrns[i] << " ";
            }
            cout << endl;
        }
    }
    // if(bb && nIter%ld != 0){
    //     Tmodp.start();
    //     Frns.reduce_modp_rnsmajor_scal_quad(rowdim*blockSize,  FFPACK::rns_double_elt_ptr(Yrns, 1));
    //     Tmodp.stop();
    //     modpTime += Tmodp.usertime();
    // }else if(!bb && nIter%ld != 0){
    //     Tmodp.start();
    //     Frns.reduce_modp_rnsmajor_scal_quad(rowdim*blockSize,  FFPACK::rns_double_elt_ptr(Xrns, 1));
    //     Tmodp.stop();
    //     modpTime += Tmodp.usertime();
    // }

    // Reconstruct Y from Yrns
    RNS.convert(rowdim*blockSize, Y.data(), Yrns);
    Ttotal.stop();
    for(size_t i = 0 ; i < rowdim*blockSize ; ++i){
        if(Y[i] < 0){
            Integer q = -Y[i] / p;
            Y[i] = p - (-Y[i] - p*q);
        }
        Y[i] %= p;
    }
    cout << "Y res:" << endl;
    for(size_t i = 0, end = (Y.size()>20)?20:Y.size() ; i < end ; ++i){
        cout << Y[i] << " ";
    }
    cout << endl;
    cout << nIter << " iterations in " << Ttotal << endl;
    cout << "spmm: " << spmmTime << endl;
    cout << "modp: " << modpTime << endl;

    fflas_delete(Xrns);
    fflas_delete(Yrns);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
