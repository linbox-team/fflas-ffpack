/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//#include "goto-def.h"

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

#include "givaro/modular.h"
#include "givaro/modular-balanced.h"

#include "fflas-ffpack/config-blas.h"
// #include "fflas-ffpac/field/modular-double.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
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

namespace details_spmv {
template <class Field> struct Coo {
  private:
    using Self = Coo<Field>;

  public:
    typename Field::Element val = 0;
    index_t col = 0;
    index_t row = 0;

    Coo() = default;
    Coo(typename Field::Element v, index_t r, index_t c) : val(v), col(c), row(r) {}
    Coo(const Self &) = default;
    Coo(Self &&) = default;

    Self &operator=(const Self &) = default;
    Self &operator=(Self &&) = default;
};
}

template <class Field>
void readSmsFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
                   typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz) {
    using namespace details_spmv;
    std::ifstream file(path, std::ios::in);
    std::vector<std::string> tokens;
    std::string line;
    // while(std::getline(file, line) && line.size()!=0);
    std::getline(file, line);
    std::istringstream is(line);
    // cout << "line : " << line << endl;
    std::copy(std::istream_iterator<std::string>(is), std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string>>(tokens));
    // cout << tokens.size() << endl;
    // cout << " " << std::stoull(tokens[0]) << " " << std::stoull(tokens[1]) << endl;
    rowdim = static_cast<index_t>(std::stoull(tokens[0]));
    coldim = static_cast<index_t>(std::stoull(tokens[1]));
    std::vector<Coo<Field>> data;
    nnz = 0;
    while (std::getline(file, line)) {
        tokens.resize(0);
        std::istringstream iss(line);

        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string>>(tokens));

        if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
            typename Field::Element v;
            f.init(v, std::stol(tokens[2]));
            index_t r = (index_t)(std::stoull(tokens[0])) - 1;
            index_t c = (index_t)(std::stoull(tokens[1])) - 1;
            data.emplace_back(v, r, c);
        }
    }
    std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });
    auto rowmax = (std::max_element(data.begin(), data.end(),
                                    [](const Coo<Field> &a, const Coo<Field> &b) { return a.row < b.row; }))->row;
    if (rowdim != rowmax + 1) {
        cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax << endl;
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

template <class Field>
void readSprFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
                   typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz) {
    using namespace details_spmv;
    std::ifstream file(path, std::ios::in);
    std::vector<std::string> tokens;
    std::string line;
    // while(std::getline(file, line) && line.size()!=0);
    std::getline(file, line);
    std::istringstream is(line);
    // cout << "line : " << line << endl;
    std::copy(std::istream_iterator<std::string>(is), std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string>>(tokens));
    // cout << tokens.size() << endl;
    // cout << " " << std::stoull(tokens[0]) << " " << std::stoull(tokens[1]) << endl;
    rowdim = static_cast<index_t>(std::stoull(tokens[0]));
    coldim = static_cast<index_t>(std::stoull(tokens[1]));
    std::vector<Coo<Field>> data;
    nnz = 0;
    uint64_t itLine = 0;
    while (std::getline(file, line)) {
        tokens.resize(0);
        std::istringstream iss(line);

        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string>>(tokens));

        // if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
        uint64_t nElements = stoull(tokens[0]);
        for (uint64_t i = 0; i < nElements; ++i) {
            index_t c = std::stoull(tokens[2 * i + 1]) - 1;
            typename Field::Element v;
            int64_t vtmp = std::stoll(tokens[2 * (i + 1)]);
            f.init(v, vtmp);
            data.emplace_back(v, itLine, c);
        }
        // typename Field::Element v;
        // f.init(v, std::stol(tokens[2]));
        // index_t r = (index_t)(std::stoull(tokens[0])) - 1;
        // index_t c = (index_t)(std::stoull(tokens[1])) - 1;
        // data.emplace_back(v, r, c);
        // }
        ++itLine;
    }
    std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });
    auto rowmax = (std::max_element(data.begin(), data.end(),
                                    [](const Coo<Field> &a, const Coo<Field> &b) { return a.row < b.row; }))->row;
    if (rowdim != rowmax + 1) {
        cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax << endl;
        rowdim = rowmax;
    }
    row = fflas_new<index_t>(data.size());
    col = fflas_new<index_t>(data.size());
    val = fflas_new(f, data.size(), 1);
    nnz = data.size();
    cout << "nnz : " << nnz << endl;
    for (size_t i = 0, end = data.size(); i < end; ++i) {
        val[i] = data[i].val;
        col[i] = data[i].col;
        row[i] = data[i].row;
    }
}

template <class Field>
void readSmsFormat2(const std::string &path, const Field &f, index_t *&row, index_t *&col,
                    typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz) {
    using namespace details_spmv;
    FILE *file;

    if ((file = fopen(path.c_str(), "r")) == NULL) {
        cout << "Error: " << path << " does not name a valid/readable file." << endl;
        return;
    }
    char code[1];
    unsigned long long int rowdim_, coldim_;
    int err;
    err = fscanf(file, "%llu %llu %s", &rowdim_, &coldim_, code);
    cout << rowdim_ << " " << coldim_ << " " << code << endl;

    rowdim = rowdim_;
    coldim = coldim_;
    std::vector<Coo<Field>> data;
    nnz = 0;

    unsigned long long int rowtmp, coltmp;
    long long int vtmp;
    while (fscanf(file, "%llu %llu %lld", &rowtmp, &coltmp, &vtmp) > 0) {
        if (vtmp != 0) {
            data.emplace_back(vtmp, rowtmp, coltmp);
            // cout << rowtmp << " " << coltmp << " " << vtmp << endl;
        }
    }

    std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });
    auto rowmax = (std::max_element(data.begin(), data.end(),
                                    [](const Coo<Field> &a, const Coo<Field> &b) { return a.row < b.row; }))->row;
    // check if empty rows at the end
    if (rowdim != rowmax + 1) {
        cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax << endl;
        rowdim = rowmax;
    }
    row = fflas_new<index_t>(data.size());
    col = fflas_new<index_t>(data.size());
    val = fflas_new(f, data.size(), 1);
    nnz = data.size();
    cout << "nnz : " << nnz << endl;
    for (size_t i = 0, end = data.size(); i < end; ++i) {
        val[i] = data[i].val;
        col[i] = data[i].col;
        row[i] = data[i].row;
    }
}

template <class MatT, class Field, class IndexT>
std::pair<double, uint64_t> test_fspmv(size_t iter, const Field &F, IndexT *row, IndexT *col,
                                       typename Field::Element_ptr dat, index_t rowdim, index_t coldim, uint64_t nnz,
                                       typename Field::Element_ptr x, typename Field::Element_ptr y,
                                       typename Field::Element beta) {
    MatT matrix;
    sparse_init(F, matrix, row, col, dat, rowdim, coldim, nnz);
    TTimer time;
    time.clear();
    time.start();
    for (size_t i = 0; i < iter; ++i)
        fspmv(F, matrix, x, 1, y);
    time.stop();
    sparse_delete(matrix);
    return make_pair(time.usertime(), matrix.nElements);
}

template <class T1, class T2, class T> void print_res(pair<T1, T2> &p, size_t iter, T as, int blocksize = 1) {
    std::cout << "Time: " << p.first / double(iter)
              << " Gflops: " << (2 * blocksize * p.second) / 1000000000. / p.first * double(iter);
    FFLAS::writeCommandString(std::cout, as) << std::endl;
}

int main(int argc, char **argv) {

    using Field = Givaro::Modular<double>;
    using Element = typename Field::Element;

    size_t iter = 10;
    int q = 1009;
    int s = 0;
    std::string matrixFile = "";

    Argument as[] = { { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INT, &q },
                      { 'i', "-i R", "Set number of repetitions.", TYPE_INT, &iter },
                      { 's', "-s S", "Compute and print matrix statistics.", TYPE_INT, &s },
                      { 'f', "-f FILE", "Set matrix file.", TYPE_STR, &matrixFile },
                      END_OF_ARGUMENTS };

    // matrixFile = "matrix/cis.mk8-8.sms";
    // matrixFile = "matrix/GL7d17.sms";
    // matrixFile = "data/mat11.sms";

    FFLAS::parseArguments(argc, argv, as);

    // cout << matrixFile << endl;

    Field F(q);

    index_t *row = nullptr, *col = nullptr;
    typename Field::Element_ptr dat;
    index_t rowdim = 0, coldim = 0;
    uint64_t nnz;

    if (matrixFile.find(".sms") != std::string::npos) {
        readSmsFormat(matrixFile, F, row, col, dat, rowdim, coldim, nnz);
    } else if (matrixFile.find(".spr") != std::string::npos) {
        readSprFormat(matrixFile, F, row, col, dat, rowdim, coldim, nnz);
    }

    if (s) {
        //auto stats = sparse_details::getStat(F, row, col, dat, rowdim, coldim, nnz);
        //std::cout << "Sparse Matrix statistics : " << std::endl;
        //stats.print();
        //std::cout << std::endl;
    }

    auto x = FFLAS::fflas_new(F, coldim, 1, Alignment::CACHE_LINE);
    auto y = FFLAS::fflas_new(F, rowdim, 1, Alignment::CACHE_LINE);

    for (size_t i = 0; i < coldim; ++i) {
        x[i] = 1;
    }

    for (size_t i = 0; i < rowdim; ++i) {
        y[i] = 0;
    }

    auto coo =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::COO>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "COO : ";
    print_res(coo, iter, as);
    auto coozo =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::COO_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "COO_ZO : ";
    print_res(coozo, iter, as);
    auto csr =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::CSR>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "CSR : ";
    print_res(csr, iter, as);
    auto csrzo =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::CSR_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "CSR_ZO : ";
    print_res(csrzo, iter, as);
    auto ell =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::ELL>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "ELL : ";
    print_res(ell, iter, as);
    auto ellzo =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::ELL_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "ELL_ZO : ";
    print_res(ellzo, iter, as);
    // auto ellsimd = test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::ELL_simd>>(iter, F, row, col, dat, rowdim, coldim,
                                                                              // nnz, x, y, 1);
    // cout << "ELL_simd : ";
    // print_res(ellsimd, iter, as);
    // auto ellsimdzo = test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::ELL_simd_ZO>>(iter, F, row, col, dat, rowdim,
    //                                                                                coldim, nnz, x, y, 1);
    // cout << "ELL_simd_ZO : ";
    // print_res(ellsimdzo, iter, as);
    auto csrhyb =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::CSR_HYB>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "CSR_HYB : ";
    print_res(csrhyb, iter, as);
    auto hybzo =
        test_fspmv<Sparse<Field, FFLAS::SparseMatrix_t::HYB_ZO>>(iter, F, row, col, dat, rowdim, coldim, nnz, x, y, 1);
    cout << "HYB_ZO : ";
    print_res(hybzo, iter, as);
    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14

    // std::cout << "Time: " << coo.first / double(iter)
    // 		  << " Gflops: " << (2*coo.second)/1000000000. / coo.first * double(iter);
    // FFLAS::writeCommandString(std::cout, as) << std::endl;

    //  std::cout << "Time: " << csr.first / double(iter)
    //        << " Gflops: " << (2*csr.second)/1000000000. / csr.first * double(iter);
    //  FFLAS::writeCommandString(std::cout, as) << std::endl;
    return 0;
}
