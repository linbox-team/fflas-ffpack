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

namespace FFLAS{

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
        std::vector<uint64_t> denseRows;
        std::vector<uint64_t> denseCols;
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
        std::vector<int64_t> rows(rowdim+1);
        std::vector<int64_t> cols(coldim);
        std::fill(rows.begin(), rows.end(), 0);
        std::fill(cols.begin(), cols.end(), 0);
        for (uint64_t i = 0; i < nnz; ++i) {
            cols[col[i]]++;
            if (F.isOne(val[i])) {
                stats.nOnes++;
            } else if (F.isMOne(val[i])) {
                stats.nMOnes++;
            } else {
                stats.nOthers++;
            }
        }
        rows[0] = row[0];
        for(size_t i = 1 ; i < rowdim+1 ; ++i){
            rows[i] = row[i] - row[i-1];
        }
        stats.nEmptyRows = std::count(rows.begin(), rows.end(), 0);
        stats.nEmptyCols = std::count(cols.begin(), cols.end(), 0);
        auto rowMinMax = std::minmax_element(rows.begin(), rows.end());
        auto colMinMax = std::minmax_element(cols.begin(), cols.end());
        stats.minRow = (*(rowMinMax.first));
        stats.maxRow = (*(rowMinMax.second));
        stats.minCol = (*(colMinMax.first));
        stats.maxCol = (*(colMinMax.second));
        stats.averageRow = std::accumulate(rows.begin(), rows.end(), 0) / rowdim;
        stats.averageCol = std::accumulate(cols.begin(), cols.end(), 0) / coldim;
        stats.deviationRow = (uint64_t)computeDeviation(rows.begin(), rows.end());
        stats.deviationCol = (uint64_t)computeDeviation(cols.begin(), cols.end());
        stats.nDenseRows = std::count_if(rows.begin(), rows.begin(),
                                         [rowdim](uint64_t &x) { return x >= DENSE_THRESHOLD * rowdim; });
        stats.nDenseCols = std::count_if(cols.begin(), cols.begin(),
                                         [coldim](uint64_t &x) { return x >= DENSE_THRESHOLD * coldim; });
        return stats;
    }

}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
