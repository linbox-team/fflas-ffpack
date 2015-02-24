/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla <bastien.vialla@lirmm.fr>
 *
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

/** @file fflas/fflas_sparse.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_utils_H
#define __FFLASFFPACK_fflas_fflas_sparse_utils_H

#include <algorithm>
#include <numeric>
#include <vector>

struct StatsMatrix {
    uint64_t rowdim = 0;
    uint64_t coldim = 0;
    uint64_t nOnes = 0;
    uint64_t nMOnes = 0;
    uint64_t nOthers = 0;
    uint64_t nnz = 0;
    uint64_t maxRow = 0;
    uint64_t minRow = 0;
    uint64_t averageRow = 0;
    uint64_t deviationRow = 0;
    uint64_t maxCol = 0;
    uint64_t minCol = 0;
    uint64_t averageCol = 0;
    uint64_t deviationCol = 0;
    uint64_t minColDifference = 0;
    uint64_t maxColDifference = 0;
    uint64_t averageColDifference = 0;
    uint64_t deviationColDifference = 0;
    uint64_t minRowDifference = 0;
    uint64_t maxRowDifference = 0;
    uint64_t averageRowDifference = 0;
    uint64_t deviationRowDifference = 0;
    uint64_t nDenseRows = 0;
    uint64_t nDenseCols = 0;
    uint64_t nEmptyRows = 0;
    uint64_t nEmptyCols = 0;
    uint64_t nEmptyColsEnd = 0;
    std::vector<uint64_t> rows;
    std::vector<uint64_t> denseRows;
    std::vector<uint64_t> denseCols;
    std::vector<uint64_t> cols;
};

template <class It> double computeDeviation(It begin, It end) {
    using T = typename std::decay<decltype(*begin)>::type;
    T average = 0;
    average = std::accumulate(begin, end, 0) / (end - begin);
    T sum = 0;
    for (It i = begin; i != end; ++i) {
        sum += ((*(i)) - average) * ((*(i)) - average);
    }
    return std::sqrt(sum / (end - begin));
}

template <class Field>
StatsMatrix getStat(const Field &F, const index_t *row, const index_t *col, typename Field::ConstElement_ptr val,
              uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
    StatsMatrix stats;
    stats.nnz = nnz;
    stats.rowdim = rowdim;
    stats.coldim = coldim;
    stats.rows.resize(rowdim);
    stats.cols.resize(coldim);
    std::fill(stats.rows.begin(), stats.rows.end(), 0);
    std::fill(stats.cols.begin(), stats.cols.end(), 0);
    for (uint64_t i = 0; i < nnz; ++i) {
        stats.rows[row[i]]++;
        stats.cols[col[i]]++;
        if (F.isOne(val[i])) {
            stats.nOnes++;
        } else if (F.isMOne(val[i])) {
            stats.nMOnes++;
        } else {
            stats.nOthers++;
        }
    }
    stats.nEmptyRows = std::count(stats.rows.begin(), stats.rows.end(), 0);
    stats.nEmptyCols = std::count(stats.cols.begin(), stats.cols.end(), 0);
    auto rowMinMax = std::minmax_element(stats.rows.begin(), stats.rows.end());
    auto colMinMax = std::minmax_element(stats.cols.begin(), stats.cols.end());
    stats.minRow = (*(rowMinMax.first));
    stats.maxRow = (*(rowMinMax.second));
    stats.minCol = (*(colMinMax.first));
    stats.maxCol = (*(colMinMax.second));
    std::vector<uint64_t> rowDiff(nnz);
    std::vector<uint64_t> colDiff(nnz);
    // std::adjacent_difference(row, row+nnz, rowDiff.begin());
    // std::adjacent_difference(col, col+nnz, colDiff.begin());
    // auto rowDiffMinMax = std::minmax_element(rowDiff.begin(), rowDiff.end());
    // auto colDiffMinMax = std::minmax_element(colDiff.begin(), colDiff.end());
    // stats.minRowDifference = *(rowDiffMinMax.first);
    // stats.maxRowDifference = *(rowDiffMinMax.second);
    // stats.minColDifference = *(colDiffMinMax.first);
    // stats.maxColDifference = *(colDiffMinMax.second);
    stats.averageRow = std::accumulate(stats.rows.begin(), stats.rows.end(), 0) / rowdim;
    stats.averageCol = std::accumulate(stats.cols.begin(), stats.cols.end(), 0) / coldim;
    // stats.averageRowDifference = std::accumulate(rowDiff.begin(), rowDiff.end(), 0)/nnz;
    // stats.averageColDifference = std::accumulate(colDiff.begin(), colDiff.end(), 0)/nnz;
    stats.deviationRow = (uint64_t)computeDeviation(stats.rows.begin(), stats.rows.end());
    stats.deviationCol = (uint64_t)computeDeviation(stats.cols.begin(), stats.cols.end());
    // stats.deviationRowDifference = (uint64_t)computeDeviation(rowDiff.begin(), rowDiff.end(),
    // stats.averageRowDifference);
    // stats.deviationColDifference = (uint64_t)computeDeviation(colDiff.begin(), colDiff.end(),
    // stats.averageColDifference);
    stats.nDenseRows = std::count_if(stats.rows.begin(), stats.rows.begin(),
                                     [rowdim](uint64_t &x) { return x >= DENSE_THRESHOLD * rowdim; });
    stats.nDenseCols = std::count_if(stats.cols.begin(), stats.cols.begin(),
                                     [coldim](uint64_t &x) { return x >= DENSE_THRESHOLD * coldim; });
    return stats;
}

#endif
