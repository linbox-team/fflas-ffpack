/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :     Bastien Vialla <bastien.vialla@lirmm.fr>
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 * ========LICENCE========
 *.
 */

#include "fflas-ffpack/fflas/fflas.h"
// #include "fflas-ffpack/fflas/fflas_fspmv.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/field/modular-double.h"
#include "fflas-ffpack/field/unparametric.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iterator>
#include <cstdlib>
#include <sstream>

using namespace FFLAS;
using namespace FFPACK;
using namespace std;

template <typename T> T from_string(std::string const & s) {
    std::stringstream ss(s);
    T result;
    ss >> result;    // TODO handle errors
    return result;
}

namespace details_spmv {
template <class Field> struct Coo {
  private:
    using Self = Coo<Field>;

  public:
    typename Field::Element val = 0;
    index_t col = 0;
    index_t row = 0;

    Coo() = default;
    Coo(typename Field::Element v, index_t r, index_t c)
        : val(v), col(c), row(r) {}
    Coo(const Self &) = default;
    Coo(Self &&) = default;

    Self &operator=(const Self &) = default;
    Self &operator=(Self &&) = default;
};
}

// TODO : faster version using fscanf
template <class Field>
void readSmsFormat(const std::string &path, const Field &f, index_t *&row,
                   index_t *&col, typename Field::Element_ptr &val,
                   index_t &rowdim, index_t &coldim, uint64_t &nnz) {
    using namespace details_spmv;
    std::ifstream file(path, std::ios::in);
    std::vector<std::string> tokens;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::copy(std::istream_iterator<std::string>(is),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string>>(tokens));
    rowdim = from_string<index_t>(tokens[0]);
    //static_cast<index_t>(stoul(tokens[0].c_str(),0,10));
    coldim = from_string<index_t>(tokens[0]);
    // static_cast<index_t>(stoul(tokens[1].c_str(),0,10));
    std::vector<Coo<Field>> data;
    nnz = 0;
    while (std::getline(file, line)) {
        tokens.resize(0);
        std::istringstream iss(line);

        std::copy(std::istream_iterator<std::string>(iss),
                  std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string>>(tokens));

        if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
            typename Field::Element v;
            long int vtmp = from_string<int64_t>(tokens[2]);
            f.init(v, vtmp);
            index_t r = from_string<index_t>(tokens[0]) - 1;
            // (index_t)(stoul(tokens[0].c_str(),0,10)) - 1;
            index_t c = from_string<index_t>(tokens[0]) - 1;
            // (index_t)(stoul(tokens[1].c_str(),0,10)) - 1;
            data.emplace_back(v, r, c);
        }
    }
    std::sort(data.begin(), data.end(),
              [](const Coo<Field> &a, const Coo<Field> &b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });
    auto rowmax =
        (std::max_element(data.begin(), data.end(),
                          [](const Coo<Field> &a, const Coo<Field> &b) {
             return a.row < b.row;
         }))->row;
    if (rowdim != rowmax + 1) {
        cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax
             << endl;
        rowdim = rowmax;
    }
    row = fflas_new<index_t>(data.size());
    col = fflas_new<index_t>(data.size());
    val = fflas_new(f, data.size(), 1);
    nnz = data.size();
    // cout << "nnz : " << nnz << endl;
    for (size_t i = 0, end = data.size(); i < end; ++i) {
        val[i] = data[i].val;
        col[i] = data[i].col;
        row[i] = data[i].row;
    }
}

template <class PtrT> void testEq(PtrT y1, PtrT y2, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
        if (y1[i] != y2[i]) {
            cout << "Error " << i << endl;
            cout << y1[i] << " != " << y2[i] << endl;
            break;
        }
    }
}

template <class MatT, class Field, class IndexT>
void test_spmv(const Field &F, IndexT *row, IndexT *col,
               typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
               uint64_t nnz, typename Field::Element_ptr x,
               typename Field::Element_ptr y, typename Field::Element beta) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    fspmv(F, matrix, x, 1, y);
    sparse_delete(matrix);
}

template <class Field, class IndexT>
void
test_spmv_sell(const Field &F, IndexT *row, IndexT *col,
               typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
               uint64_t nnz, Sparse<Field, SparseMatrix_t::SELL> &matrix,
               typename Field::Element_ptr x, typename Field::Element_ptr y,
               typename Field::Element beta) {
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    fspmv(F, matrix, x, 1, y);
    auto tmp = fflas_new(F, rowdim, 1);
    for (size_t i = 0; i < rowdim; ++i) {
        tmp[i] = y[matrix.perm[i]];
    }
    for (size_t i = 0; i < rowdim; ++i) {
        y[i] = tmp[i];
    }
    sparse_delete(matrix);
    fflas_delete(tmp);
}

template <class MatT, class Field, class IndexT>
void test_spmm(const Field &F, IndexT *row, IndexT *col,
               typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
               uint64_t nnz, int blockSize, typename Field::Element_ptr x,
               int ldx, typename Field::Element_ptr y, int ldy,
               typename Field::Element beta) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    fspmm(F, matrix, blockSize, x, ldx, beta, y, ldy);
    sparse_delete(matrix);
}

template <class MatT, class Field, class IndexT>
void test_pspmm(const Field &F, IndexT *row, IndexT *col,
                typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
                uint64_t nnz, int blockSize, typename Field::Element_ptr x,
                int ldx, typename Field::Element_ptr y, int ldy,
                typename Field::Element beta) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    pfspmm(F, matrix, blockSize, x, ldx, beta, y, ldy);
    sparse_delete(matrix);
}

template <class MatT, class Field, class IndexT>
void test_pspmv(const Field &F, IndexT *row, IndexT *col,
                typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
                uint64_t nnz, typename Field::Element_ptr x,
                typename Field::Element_ptr y, typename Field::Element beta) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    pfspmv(F, matrix, x, 1, y);
    sparse_delete(matrix);
}

template <class Field, class IndexT>
void
test_pspmv_sell(const Field &F, IndexT *row, IndexT *col,
                typename Field::Element_ptr dat, index_t rowdim, index_t coldim,
                uint64_t nnz, Sparse<Field, SparseMatrix_t::SELL> &matrix,
                typename Field::Element_ptr x, typename Field::Element_ptr y,
                typename Field::Element beta) {
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    pfspmv(F, matrix, x, 1, y);
    auto tmp = fflas_new(F, rowdim, 1);
    for (size_t i = 0; i < rowdim; ++i) {
        tmp[i] = y[matrix.perm[i]];
    }
    for (size_t i = 0; i < rowdim; ++i) {
        y[i] = tmp[i];
    }
    sparse_delete(matrix);
    fflas_delete(tmp);
}

int main(int argc, char **argv) {
    // using Field = UnparametricField<double>;
    using Field = Modular<float>;
    // using Element = typename Field::Element;

    Field F(101);
    int nbTests = 25;
    std::string path;

    index_t *row = nullptr, *col = nullptr;
    typename Field::Element_ptr dat;

    index_t rowdim, coldim;
    uint64_t nnz;

    // if(argc > 1)
    //     path = argv[1];

    path = "data/mat11.sms";

    readSmsFormat(path, F, row, col, dat, rowdim, coldim, nnz);

    auto x = fflas_new(F, coldim, 1, Alignment::CACHE_LINE);
    auto y = fflas_new(F, rowdim, 1, Alignment::CACHE_LINE);
    auto y1 = fflas_new(F, rowdim, 1, Alignment::CACHE_LINE);

    for (size_t i = 0; i < coldim; ++i) {
        x[i] = 1;
    }

    for (size_t i = 0; i < rowdim; ++i) {
        y[i] = 0;
        y1[i] = 0;
    }

    /************************************************************************************
     *
     * SPMV
     *
     *************************************************************************************/
    cout << "=== spmv ===" << endl;

    test_spmv<Sparse<Field, SparseMatrix_t::CSR>>(F, row, col, dat, rowdim,
                                                  coldim, nnz, x, y, 1);
    cout << "CSR: OK" << endl;

    test_spmv<Sparse<Field, SparseMatrix_t::COO>>(F, row, col, dat, rowdim,
                                                  coldim, nnz, x, y1, 1);

    cout << "COO: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_spmv<Sparse<Field, SparseMatrix_t::ELL>>(F, row, col, dat, rowdim,
                                                  coldim, nnz, x, y1, 1);

    cout << "ELL: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_spmv<Sparse<Field, SparseMatrix_t::ELL_simd>>(F, row, col, dat, rowdim,
                                                       coldim, nnz, x, y1, 1);

    cout << "ELL_simd: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_spmv<Sparse<Field, SparseMatrix_t::CSR_HYB>>(F, row, col, dat, rowdim,
                                                      coldim, nnz, x, y1, 1);

    cout << "CSR_HYB: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_spmv<Sparse<Field, SparseMatrix_t::HYB_ZO>>(F, row, col, dat, rowdim,
                                                     coldim, nnz, x, y1, 1);

    cout << "HYB_ZO: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    Sparse<Field, SparseMatrix_t::SELL> A;
    test_spmv_sell(F, row, col, dat, rowdim, coldim, nnz, A, x, y1, 1);

    cout << "SELL: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    /************************************************************************************
     *
     * pSPMV
     *
     *************************************************************************************/

    cout << "=== pspmv ===" << endl;

    test_pspmv<Sparse<Field, SparseMatrix_t::CSR>>(F, row, col, dat, rowdim,
                                                   coldim, nnz, x, y1, 1);
    cout << "CSR: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_pspmv<Sparse<Field, SparseMatrix_t::ELL>>(F, row, col, dat, rowdim,
                                                   coldim, nnz, x, y1, 1);

    cout << "ELL: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_pspmv<Sparse<Field, SparseMatrix_t::ELL_simd>>(
        F, row, col, dat, rowdim, coldim, nnz, x, y1, 1);

    cout << "ELL_simd: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_pspmv<Sparse<Field, SparseMatrix_t::CSR_HYB>>(F, row, col, dat, rowdim,
                                                       coldim, nnz, x, y1, 1);

    cout << "CSR_HYB: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    test_pspmv<Sparse<Field, SparseMatrix_t::HYB_ZO>>(F, row, col, dat, rowdim,
                                                      coldim, nnz, x, y1, 1);

    cout << "HYB_ZO: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    Sparse<Field, SparseMatrix_t::SELL> A1;
    test_pspmv_sell(F, row, col, dat, rowdim, coldim, nnz, A1, x, y1, 1);

    cout << "SELL: " << ((std::equal(y, y + rowdim, y1)) ? "OK" : "ERROR")
         << endl;

    for (size_t i = 0; i < rowdim; ++i) {
        y1[i] = 0;
    }

    // // test_spmm<Sparse<Field, SparseMatrix_t::CSR>>(F, row, col, dat,
    // rowdim,
    // coldim, nnz, 1, x, 1, y, 1, 1);
    // // test_pspmm<Sparse<Field, SparseMatrix_t::CSR>>(F, row, col, dat,
    // rowdim,
    // coldim, nnz, 1, x, 1, y, 1, 1);
    // // test_spmm<Sparse<Field, SparseMatrix_t::COO_ZO>>(F, row, col, dat,
    // rowdim, coldim, nnz, 1, x, 1, y, 1, 1);
    // // test_spmm<Sparse<Field, SparseMatrix_t::CSR>>(F, row, col, dat,
    // rowdim,
    // coldim, nnz, 1, x, 1, y, 1, 1);
    // // test_spmv<Sparse<Field, SparseMatrix_t::ELL_ZO>>(F, row, col, dat,
    // rowdim, coldim, nnz, x, y1, 1);

    // for(size_t i = 0 ; i < 11 ; ++i)
    // {
    //     cout << y[i] << " ";
    // }
    // cout << endl;

    // for(size_t i = 0 ; i < 11 ; ++i)
    // {
    //     cout << y1[i] << " ";
    // }
    // cout << endl;

    // auto bb = std::equal(y, y+rowdim, y1);

    // cout << ((bb) ? "CORRECT" : "ERROR") << endl;

    // if(!bb)
    //     testEq(y, y1, rowdim);

    fflas_delete(x);
    fflas_delete(y);
    fflas_delete(y1);

    return 0;
}
