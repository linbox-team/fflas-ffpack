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

#include "fflas-ffpack/config-blas.h"
// #include "fflas-ffpac/field/modular-double.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_OPENMP
typedef FFLAS::OMPTimer TTimer;
#else
typedef FFLAS::Timer TTimer;
#endif

using namespace std;
using namespace FFLAS;

template <typename T> T from_string(std::string const &s) {
    std::stringstream ss(s);
    T result;
    ss >> result; // TODO handle errors
    return result;
}


template <class MatT, class Field, class IndexT>
std::pair<double, uint64_t> test_fspmm(size_t iter, const Field &F, IndexT *row, IndexT *col,
                                       typename Field::Element_ptr dat, index_t rowdim, index_t coldim, uint64_t nnz,
                                       int blocksize, typename Field::Element_ptr x, int ldx,
                                       typename Field::Element beta, typename Field::Element_ptr y, int ldy) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    TTimer time;
    time.clear();
    time.start();
    for (size_t i = 0; i < iter; ++i)
        fspmm(F, matrix, blocksize, x, ldx, 1, y, ldy);
    time.stop();
    sparse_delete(matrix);
    return make_pair(time.usertime(), matrix.nElements);
}

template <class T1, class T2, class T> void print_res(pair<T1, T2> &p, size_t iter, T as, int blocksize) {
    //	cout << 2*p.second*blocksize*iter << endl;
    std::cout << "Time: " << p.first / double(iter)
    << " Gfops: " << ((2*blocksize*p.second)/(1000000.*p.first))*(double(iter)/1000) ;
    FFLAS::writeCommandString(std::cout, as) << std::endl;
}

int main(int argc, char **argv) {

    using Field = Givaro::Modular<double>;
    using Element = typename Field::Element;

    size_t iter = 10;
    int q = 1009;
    int blocksize = 4;
    int s = 0;
    std::string matrixFile = "";

    Argument as[] = { { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INT, &q },
        { 'b', "-b Q", "Set the block size.", TYPE_INT, &blocksize },
        { 'i', "-i R", "Set number of repetitions.", TYPE_INT, &iter },
        { 's', "-s S", "Compute and print matrix statistics.", TYPE_INT, &s },
        { 'f', "-f FILE", "Set matrix file.", TYPE_STR, &matrixFile },
        END_OF_ARGUMENTS };

    // matrixFile = "matrix/cis.mk8-8.sms";
    // matrixFile = "matrix/M06-D9.sms";
    // matrixFile = "matrix/GL7d17.sms";
    // matrixFile = "data/mat11.sms";

    FFLAS::parseArguments(argc, argv, as);

    // cout << matrixFile << endl;

    Field F(q);

    index_t *row = nullptr, *col = nullptr;
    typename Field::Element_ptr dat;
    index_t rowdim, coldim;
    uint64_t nnz;

    index_t * st = nullptr ;
    readSmsFormat(matrixFile, F, st, col, dat, rowdim, coldim, nnz);
    row = fflas_new<index_t>(nnz);
    for (index_t j = 0 ; j < rowdim ; ++j) {
        for (index_t k = st[j] ; k < st[j+1] ; ++k)
            row[k] = j ;
    }

    if (s) {
        // auto stats = sparse_details::getStat(F, row, col, dat, rowdim, coldim, nnz);
        // std::cout << "Sparse Matrix statistics : " << std::endl;
        // stats.print();
        std::cout << std::endl;
    }

    auto x = FFLAS::fflas_new(F, coldim, blocksize, Alignment::CACHE_LINE);
    auto y = FFLAS::fflas_new(F, rowdim, blocksize, Alignment::CACHE_LINE);

    for (size_t i = 0; i < coldim * blocksize; ++i) {
        x[i] = 1;
    }

    for (size_t i = 0; i < rowdim * blocksize; ++i) {
        y[i] = 0;
    }

    // auto coo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::COO>>(iter, F, row, col, dat, rowdim, coldim, nnz,
    // blocksize, x, blocksize, 1, y, blocksize);
    // cout << "COO : ";
    // print_res(coo, iter, as);

    // auto coozo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::COO_ZO>>(iter, F, row, col, dat, rowdim, coldim,
    // nnz, blocksize, x, blocksize, 1, y, blocksize);
    // cout << "COO_ZO : ";
    // print_res(coozo, iter, as);
    auto csr = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::CSR>>(iter, F, row, col, dat, rowdim, coldim, nnz,
                                                                     blocksize, x, blocksize, 1, y, blocksize);
    cout << "CSR : ";
    print_res(csr, iter, as, blocksize);
    auto ell = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::ELL>>(iter, F, row, col, dat, rowdim, coldim, nnz,
                                                                     blocksize, x, blocksize, 1, y, blocksize);
    cout << "ELL : ";
    print_res(ell, iter, as, blocksize);
    auto ellzo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::ELL_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz,
                                                                          blocksize, x, blocksize, 1, y, blocksize);
    cout << "ELL_ZO : ";
    print_res(ellzo, iter, as, blocksize);
    // auto csrzo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::CSR_ZO>>(iter, F, row, col, dat, rowdim, coldim,
    // nnz, blocksize, x, blocksize, 1, y, blocksize);
    // cout << "CSR_ZO : ";
    // print_res(csrzo, iter, as);
    // auto ell = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::ELL>>(iter, F, row, col, dat, rowdim, coldim, nnz,
    // blocksize, x, blocksize, 1, y, blocksize);
    // cout << "ELL : ";
    // print_res(ell, iter, as);
    // auto ellzo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::ELL_ZO>>(iter, F, row, col, dat, rowdim, coldim,
    // nnz, blocksize, x, blocksize, 1, y, blocksize);
    // cout << "ELL_ZO : ";
    // print_res(ellzo, iter, as);
    auto hybzo = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::HYB_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz,
                                                                          blocksize, x, blocksize, 1, y, blocksize);
    cout << "HYB_ZO : ";
    print_res(hybzo, iter, as, blocksize);
    auto csrhyb = test_fspmm<Sparse<Field, FFLAS::SparseMatrix_t::CSR_HYB>>(iter, F, row, col, dat, rowdim, coldim, nnz,
                                                                            blocksize, x, blocksize, 1, y, blocksize);
    cout << "CSR_HYB : ";
    print_res(csrhyb, iter, as, blocksize);
    // for (size_t i = 0; i < 10*blocksize; ++i) {
    //   std::cout << y[i] << " ";
    // }
    // std::cout << std::endl;

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14

    // std::cout << "Time: " << coo.first / double(iter)
    // 		  << " Gfops: " << (2*coo.second)/1000000000. / coo.first * double(iter);
    // FFLAS::writeCommandString(std::cout, as) << std::endl;

    //  std::cout << "Time: " << csr.first / double(iter)
    //        << " Gfops: " << (2*csr.second)/1000000000. / csr.first * double(iter);
    //  FFLAS::writeCommandString(std::cout, as) << std::endl;
    fflas_delete(x);
    fflas_delete(y);
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
