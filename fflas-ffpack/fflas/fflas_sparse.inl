/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
 *              Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_sparse.inl
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_INL
#define __FFLASFFPACK_fflas_fflas_sparse_INL

namespace FFLAS {
namespace sparse_details {
template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::ModularTag) {
    if (b != 1) {
        if (b == 0) {
            for (size_t i = 0; i < m; ++i)
                y[i] = 0;
        } else if (b == -1) {
            for (size_t i = 0; i < m; ++i)
                y[i] *= -1;
        } else {
            fscalin(F, m, b, y, 1);
        }
    }
}

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::ModularTag) {
    if (b != 1) {
        if (b == 0) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    y[i * ldy + j] = 0;
                }
            }
        } else if (b == -1) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    y[i * ldy + j] *= -1;
                }
            }
        } else {
            fscalin(F, m, n, y, ldy);
        }
    }
}

template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::UnparametricTag) {
    if (b != 1) {
        if (b == 0) {
            for (size_t i = 0; i < m; ++i)
                y[i] = 0;
        } else if (b == -1) {
            for (size_t i = 0; i < m; ++i)
                y[i] *= -1;
        } else {
            fscalin(F, m, b, y, 1);
        }
    }
}

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::UnparametricTag) {
    if (b != 1) {
        if (b == 0) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    y[i * ldy + j] = 0;
                }
            }
        } else if (b == -1) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    y[i * ldy + j] *= -1;
                }
            }
        } else {
            fscalin(F, m, n, y, ldy);
        }
    }
}

template <class Field>
inline void init_y(const Field &F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y,
                   FieldCategories::GenericTag) {
    if (!F.isOne(b)) {
        if (F.isZero(b)) {
            for (size_t i = 0; i < m; ++i)
                F.assign(y[i], F.zero);
        } else if (F.isMOne(b)) {
            for (size_t i = 0; i < m; ++i)
                F.negin(y[i]);
        } else {
            fscalin(F, m, b, y, 1);
        }
    }
}

template <class Field>
inline void init_y(const Field &F, const size_t m, const size_t n, const typename Field::Element b,
                   typename Field::Element_ptr y, const int ldy, FieldCategories::GenericTag) {
    if (b != 1) {
        if (b == 0) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    F.assign(y[i * ldy + j], F.zero);
                }
            }
        } else if (b == -1) {
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    F.negin(y[i * ldy + j]);
                }
            }
        } else {
            fscalin(F, m, n, y, ldy);
        }
    }
}
} // sparse_details

namespace sparse_details {
// template <class It> double computeDeviation(It begin, It end) {
//     using T = typename std::decay<decltype(*begin)>::type;
//     T average = 0;
//     average = std::accumulate(begin, end, 0) / (end - begin);
//     T sum = 0;
//     for (It i = begin; i != end; ++i) {
//         sum += ((*(i)) - average) * ((*(i)) - average);
//     }
//     return std::sqrt(sum / (end - begin));
// }

// template <class Field>
// Stats getStat(const Field &F, const index_t *row, const index_t *col, typename Field::ConstElement_ptr val,
//               uint64_t rowdim, uint64_t coldim, uint64_t nnz) {
//     Stats stats;
//     stats.nnz = nnz;
//     stats.rowdim = rowdim;
//     stats.coldim = coldim;
//     stats.rows.resize(rowdim);
//     stats.cols.resize(coldim);
//     std::fill(stats.rows.begin(), stats.rows.end(), 0);
//     std::fill(stats.cols.begin(), stats.cols.end(), 0);
//     for (uint64_t i = 0; i < nnz; ++i) {
//         stats.rows[row[i]]++;
//         stats.cols[col[i]]++;
//         if (F.isOne(val[i])) {
//             stats.nOnes++;
//         } else if (F.isMOne(val[i])) {
//             stats.nMOnes++;
//         } else {
//             stats.nOthers++;
//         }
//     }
//     stats.nEmptyRows = std::count(stats.rows.begin(), stats.rows.end(), 0);
//     stats.nEmptyCols = std::count(stats.cols.begin(), stats.cols.end(), 0);
//     auto rowMinMax = std::minmax_element(stats.rows.begin(), stats.rows.end());
//     auto colMinMax = std::minmax_element(stats.cols.begin(), stats.cols.end());
//     stats.minRow = (*(rowMinMax.first));
//     stats.maxRow = (*(rowMinMax.second));
//     stats.minCol = (*(colMinMax.first));
//     stats.maxCol = (*(colMinMax.second));
//     std::vector<uint64_t> rowDiff(nnz);
//     std::vector<uint64_t> colDiff(nnz);
//     // std::adjacent_difference(row, row+nnz, rowDiff.begin());
//     // std::adjacent_difference(col, col+nnz, colDiff.begin());
//     // auto rowDiffMinMax = std::minmax_element(rowDiff.begin(), rowDiff.end());
//     // auto colDiffMinMax = std::minmax_element(colDiff.begin(), colDiff.end());
//     // stats.minRowDifference = *(rowDiffMinMax.first);
//     // stats.maxRowDifference = *(rowDiffMinMax.second);
//     // stats.minColDifference = *(colDiffMinMax.first);
//     // stats.maxColDifference = *(colDiffMinMax.second);
//     stats.averageRow = std::accumulate(stats.rows.begin(), stats.rows.end(), 0) / rowdim;
//     stats.averageCol = std::accumulate(stats.cols.begin(), stats.cols.end(), 0) / coldim;
//     // stats.averageRowDifference = std::accumulate(rowDiff.begin(), rowDiff.end(), 0)/nnz;
//     // stats.averageColDifference = std::accumulate(colDiff.begin(), colDiff.end(), 0)/nnz;
//     stats.deviationRow = (uint64_t)computeDeviation(stats.rows.begin(), stats.rows.end());
//     stats.deviationCol = (uint64_t)computeDeviation(stats.cols.begin(), stats.cols.end());
//     // stats.deviationRowDifference = (uint64_t)computeDeviation(rowDiff.begin(), rowDiff.end(),
//     // stats.averageRowDifference);
//     // stats.deviationColDifference = (uint64_t)computeDeviation(colDiff.begin(), colDiff.end(),
//     // stats.averageColDifference);
//     stats.nDenseRows = std::count_if(stats.rows.begin(), stats.rows.begin(),
//                                      [rowdim](uint64_t &x) { return x >= DENSE_THRESHOLD * rowdim; });
//     stats.nDenseCols = std::count_if(stats.cols.begin(), stats.cols.begin(),
//                                      [coldim](uint64_t &x) { return x >= DENSE_THRESHOLD * coldim; });
//     return stats;
// }
}

namespace sparse_details {
// non ZO matrix
template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::GenericTag, std::false_type) {
    sparse_details_impl::fspmv(F, A, x, y, FieldCategories::GenericTag());
}

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::UnparametricTag, std::false_type) {
    sparse_details_impl::fspmv(F, A, x, y, FieldCategories::UnparametricTag());
}

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::ModularTag, std::false_type) {
    if (A.delayed) {
        sparse_details::fspmv(F, A, x, y, FieldCategories::UnparametricTag(), std::false_type());
        freduce(F, A.m, y, 1);
    } else {
        sparse_details_impl::fspmv(F, A, x, y, A.kmax);
    }
}

// template <class Field, class SM>
// inline void fspmv_task(const Field &F, const index_t start, const indext_t size, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
//                   FieldCategories::GenericTag, std::false_type) {
//     sparse_details_impl::fspmv_task(F, start, size, A, x, y, FieldCategories::GenericTag());
// }

// template <class Field, class SM>
// inline void fspmv_task(const Field &F, const index_t start, const indext_t size, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
//                   FieldCategories::UnparametricTag, std::false_type) {
//     sparse_details_impl::fspmv_task(F, start, size, A, x, y, FieldCategories::UnparametricTag());
// }

// template <class Field, class SM>
// inline void fspmv_task(const Field &F, const index_t start, const indext_t size, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
//                   FieldCategories::ModularTag, std::false_type) {
//     if (A.delayed) {
//         sparse_details::fspmv_task(F, start, size, A, x, y, FieldCategories::UnparametricTag(), std::false_type());
//         freduce(F, A.m, y, 1);
//     } else {
//         sparse_details_impl::fspmv_task(F, start, size, A, x, y, A.kmax);
//     }
// }

// ZO matrix
template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::GenericTag, std::true_type) {
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::GenericTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone(F, A, x, y, FieldCategories::GenericTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::GenericTag());
        fflas_delete(x1);
    }
}

template <class Field>
inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A, typename Field::ConstElement_ptr x,
                  typename Field::Element_ptr y, FieldCategories::UnparametricTag, std::true_type) {
#ifdef __FFLASFFPACK_USE_SIMD
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#else
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#endif // SIMD
}

template <class Field>
inline void fspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_simd> &A, typename Field::ConstElement_ptr x,
                  typename Field::Element_ptr y, FieldCategories::UnparametricTag, std::true_type) {
#ifdef __FFLASFFPACK_USE_SIMD
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#else
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#endif // SIMD
}

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::UnparametricTag, std::true_type) {
    if (A.cst == 1) {
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::fspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::fspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
}

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                  FieldCategories::ModularTag, std::true_type) {
    sparse_details::fspmv<Field, SM>(F, A, x, y, FieldCategories::UnparametricTag(), std::true_type());
    freduce(F, A.m, y, 1);
}

/***************************** pfspmv ******************************/

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::GenericTag, std::false_type) {
    sparse_details_impl::pfspmv(F, A, x, y, FieldCategories::GenericTag());
}

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::UnparametricTag, std::false_type) {
    sparse_details_impl::pfspmv(F, A, x, y, FieldCategories::UnparametricTag());
}

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::ModularTag, std::false_type) {
    if (A.delayed) {
        sparse_details::pfspmv(F, A, x, y, FieldCategories::UnparametricTag(), std::false_type());
        freduce(F, A.m, y, 1);
    } else {
        sparse_details_impl::pfspmv(F, A, x, y, A.kmax);
    }
}

// ZO matrix
template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::GenericTag, std::true_type) {
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::GenericTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone(F, A, x, y, FieldCategories::GenericTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::GenericTag());
        fflas_delete(x1);
    }
}

template <class Field>
inline void pfspmv(const Field &F, const Sparse<Field, SparseMatrix_t::SELL> &A, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y, FieldCategories::UnparametricTag, std::true_type) {
#ifdef __FFLASFFPACK_USE_SIMD
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#else
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#endif // SIMD
}

template <class Field>
inline void pfspmv(const Field &F, const Sparse<Field, SparseMatrix_t::ELL_simd> &A, typename Field::ConstElement_ptr x,
                   typename Field::Element_ptr y, FieldCategories::UnparametricTag, std::true_type) {
#ifdef __FFLASFFPACK_USE_SIMD
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone_simd(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one_simd(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#else
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#endif // SIMD
}

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::UnparametricTag, std::true_type) {
    if (A.cst == 1) {
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
    } else if (A.cst == -1) {
        sparse_details_impl::pfspmv_mone(F, A, x, y, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
        fscal(F, A.n, A.cst, x, 1, x1, 1);
        sparse_details_impl::pfspmv_one(F, A, x, y, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
}

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                   FieldCategories::ModularTag, std::true_type) {
    sparse_details::pfspmv<Field, SM>(F, A, x, y, FieldCategories::UnparametricTag(), std::true_type());
    freduce(F, A.m, y, 1);
}

// /***************************** fspmm *****************************/

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, std::false_type) {
    // std::cout << "no ZO Generic" << endl;
    sparse_details_impl::fspmm(F, A, blockSize, x, ldx, y, ldy, FieldCategories::GenericTag());
}

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, std::false_type) {
// std::cout << "no ZO Unparametric" << endl;
#ifdef __FFLASFFPACK_USE_SIMD
    using simd = Simd<typename Field::Element>;
    if (((uint64_t)y % simd::alignment == 0) && ((uint64_t)x % simd::alignment == 0) &&
        (blockSize % simd::vect_size == 0)) {
        // std::cout << "no ZO Unparametric algined" << endl;
        sparse_details_impl::fspmm_simd_aligned(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
    }
// else{
//     sparse_details_impl::fspmm_simd_unaligned(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
// }
#else
    sparse_details_impl::fspmm(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
#endif
}

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, std::false_type) {
    // std::cout << "no ZO Modular" << endl;
    if (A.delayed) {
        sparse_details::fspmm(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag(),
                              typename std::false_type());
        freduce(F, A.m, blockSize, y, ldy);
    } else {
#ifdef __FFLASFFPACK_USE_SIMD
        using simd = Simd<typename Field::Element>;
        if (((uint64_t)y % simd::alignment == 0) && ((uint64_t)x % simd::alignment == 0) &&
            (blockSize % simd::vect_size == 0)) {
            sparse_details_impl::fspmm_simd_aligned(F, A, blockSize, x, ldx, y, ldy, A.kmax);
        } else {
            sparse_details_impl::fspmm_simd_unaligned(F, A, blockSize, x, ldx, y, ldy, A.kmax);
        }
#else
        sparse_details_impl::fspmm(F, A, blockSize, x, ldx, y, ldy, A.kmax);
#endif
    }
}

// ZO matrix
template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag, std::true_type) {
    // std::cout << "ZO Generic" << endl;
    if (F.isOne(A.cst)) {
        sparse_details_impl::fspmm_one(F, A, blockSize, x, ldx, y, ldy, FieldCategories::GenericTag());
    } else if (F.isMOne(A.cst)) {
        sparse_details_impl::fspmm_mone(F, A, blockSize, x, ldx, y, ldy, FieldCategories::GenericTag());
    } else {
        auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
        fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
        sparse_details_impl::fspmm_one(F, A, blockSize, x, ldx, y, ldy, FieldCategories::GenericTag());
        fflas_delete(x1);
    }
}

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag, std::true_type) {
// std::cout << "ZO Unparametric" << endl;
#ifdef __FFLASFFPACK_USE_SIMD
    using simd = Simd<typename Field::Element>;
    if (F.isOne(A.cst)) {
        if (((uint64_t)y % simd::alignment == 0) && ((uint64_t)x % simd::alignment == 0) &&
            (blockSize % simd::vect_size == 0)) {
            // std::cout << "ZO Unparametric aligned" << endl;
            sparse_details_impl::fspmm_one_simd_aligned(F, A, blockSize, x, ldx, y, ldy,
                                                        FieldCategories::UnparametricTag());
        } else {
            // std::cout << "ZO Unparametric unaligned" << endl;
            sparse_details_impl::fspmm_one_simd_unaligned(F, A, blockSize, x, ldx, y, ldy,
                                                          FieldCategories::UnparametricTag());
        }
    } else if (F.isMOne(A.cst)) {
        if (((uint64_t)y % simd::alignment == 0) && ((uint64_t)x % simd::alignment == 0) &&
            (blockSize % simd::vect_size == 0)) {
            sparse_details_impl::fspmm_mone_simd_aligned(F, A, blockSize, x, ldx, y, ldy,
                                                         FieldCategories::UnparametricTag());
        } else {
            sparse_details_impl::fspmm_mone_simd_unaligned(F, A, blockSize, x, ldx, y, ldy,
                                                           FieldCategories::UnparametricTag());
        }
    } else {
        auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
        fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
        if (((uint64_t)y % simd::alignment == 0) && ((uint64_t)x % simd::alignment == 0) &&
            (blockSize % simd::vect_size == 0)) {
            sparse_details_impl::fspmm_one_simd_aligned(F, A, blockSize, x, ldx, y, ldy,
                                                        FieldCategories::UnparametricTag());
        } else {
            sparse_details_impl::fspmm_one_simd_unaligned(F, A, blockSize, x, ldx, y, ldy,
                                                          FieldCategories::UnparametricTag());
        }
        fflas_delete(x1);
    }
#else
    if (F.isOne(A.cst)) {
        sparse_details_impl::fspmm_one(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
    } else if (F.isMOne(A.cst)) {
        sparse_details_impl::fspmm_mone(F, A, blockSize, x, ldx, y, ldy, FieldCategories::UnparametricTag());
    } else {
        auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
        fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
        sparse_details_impl::fspmm_one(F, A, blockSize, x1, ldx, y, ldy, FieldCategories::UnparametricTag());
        fflas_delete(x1);
    }
#endif
}

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  typename Field::Element_ptr y, int ldy, FieldCategories::ModularTag, std::true_type) {
    // std::cout << "ZO Modular" << endl;
    if (A.delayed) {
        sparse_details::fspmm(F, A, blockSize, x, ldx, y, ldy, typename FieldCategories::UnparametricTag(),
                              typename std::true_type());
        freduce(F, blockSize, A.m, y, ldy);
    } else {
        sparse_details_impl::fspmm(F, A, blockSize, x, ldx, y, ldy, A.kmax);
    }
}

// /***************************** fspmm *****************************/

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::GenericTag, std::false_type) {
//     if (ldx == blockSize && ldy == blockSize) {
//         sparse_details_impl::pfspmm(F, A, blockSize, x, y,
//                                     FieldCategories::GenericTag());
//     } else {
//         sparse_details_impl::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                     FieldCategories::GenericTag());
//     }
// }

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::UnparametricTag, std::false_type) {
// #ifdef __FFLASFFPACK_USE_SIMD
//     using simd = Simd<typename Field::Element>;
//     if (ldx == blockSize && ldy == blockSize) {
//         if (blockSize % simd::vect_size) {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, y, simd::loadu,
//                                         simd::storeu,
//                                         FieldCategories::UnparametricTag());
//         } else {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, y, simd::load,
//                                         simd::store,
//                                         FieldCategories::UnparametricTag());
//         }
//     } else {
//         if (blockSize % simd::vect_size) {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                         simd::loadu, simd::storeu,
//                                         FieldCategories::UnparametricTag());
//         } else {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                         simd::load, simd::store,
//                                         FieldCategories::UnparametricTag());
//         }
//     }
// #else
//     if (ldx == blockSize && ldy == blockSize) {
//         sparse_details_impl::pfspmm(F, A, blockSize, x, y,
//                                     FieldCategories::UnparametricTag());
//     } else {
//         sparse_details_impl::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                     FieldCategories::UnparametricTag());
//     }
// #endif
// }

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::ModularTag, std::false_type) {
//     if (A.delayed) {
//         sparse_details::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                FieldCategories::UnparametricTag(),
//                                typename std::false_type());
//         freduce(F, A.m, blockSize, y, ldy);
//     } else {
//         if (ldx == blockSize && ldy == blockSize) {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, y, A.kmax);
//         } else {
//             sparse_details_impl::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                         A.kmax);
//         }
//     }
// }

// // ZO matrix
// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::GenericTag, std::true_type) {
//     using Element = typename Field::Element;
//     if (ldx == blockSize && ldy == blockSize) {
//         if (F.isOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, y,
//                 [&F](Element &a, const Element &b) { F.addin(a, b); });
//         } else if (F.isMOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, y,
//                 [&F](Element &a, const Element &b) { F.subin(a, b); });
//         } else {
//             auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//             fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x1, y,
//                 [&F](Element &a, const Element &b) { F.addin(a, b); });
//             fflas_delete(x1);
//         }
//     } else {
//         if (F.isOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, ldx, y, ldy,
//                 [&F](Element &a, const Element &b) { F.addin(a, b); });
//         } else if (F.isMOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, ldx, y, ldy,
//                 [&F](Element &a, const Element &b) { F.subin(a, b); });
//         } else {
//             auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//             fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, ldx, y, ldy,
//                 [&F](Element &a, const Element &b) { F.addin(a, b); });
//             fflas_delete(x1);
//         }
//     }
// }

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::UnparametricTag, std::true_type) {
//     using Element = typename Field::Element;
// #ifdef __FFLASFFPACK_USE_SIMD
//     using simd = Simd<typename Field::Element>;
//     if (blockSize % simd::vect_size) {
//         if (ldx == blockSize && ldy == blockSize) {
//             if (F.isOne(A.cst)) {
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::add,
//                     [](Element &a, const Element &b) { a += b; }, simd::loadu,
//                     simd::storeu);
//             } else if (F.isMOne(A.cst)) {
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::sub,
//                     [](Element &a, const Element &b) { a -= b; }, simd::loadu,
//                     simd::storeu);
//             } else {
//                 auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//                 fscal(F, A.m * blockSize, A.cst, x, ldx, x1, 1);
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::add,
//                     [](Element &a, const Element &b) { a += b; }, simd::loadu,
//                     simd::storeu);
//                 fflas_delete(x1);
//             }
//         } else {
//             if (F.isOne(A.cst)) {
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::add,
//                     [](Element &a, const Element &b) { a += b; }, simd::load,
//                     simd::store);
//             } else if (F.isMOne(A.cst)) {
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::sub,
//                     [](Element &a, const Element &b) { a -= b; }, simd::load,
//                     simd::store);
//             } else {
//                 auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//                 fscal(F, A.m * blockSize, A.cst, x, ldx, x1, 1);
//                 sparse_details_impl::pfspmm_zo(
//                     F, A, blockSize, x, y, simd::add,
//                     [](Element &a, const Element &b) { a += b; }, simd::load,
//                     simd::store);
//                 fflas_delete(x1);
//             }
//         }
//     } else {
//     }
// #else
//     if (ldx == blockSize && ldy == blockSize) {
//         if (F.isOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, y,
//                 [](Element &a, const Element &b) { a += b; });
//         } else if (F.isMOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, y,
//                 [](Element &a, const Element &b) { a -= b; });
//         } else {
//             auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//             fscal(F, A.m * blockSize, A.cst, x, ldx, x1, 1);
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, y,
//                 [](Element &a, const Element &b) { a += b; });
//             fflas_delete(x1);
//         }
//     } else {
//         if (F.isOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, ldx, y, ldy,
//                 [](Element &a, const Element &b) { a += b; });
//         } else if (F.isMOne(A.cst)) {
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x, ldx, y, ldy,
//                 [](Element &a, const Element &b) { a -= b; });
//         } else {
//             auto x1 = fflas_new(F, A.m, blockSize, Alignment::CACHE_LINE);
//             fscal(F, A.m, blockSize, A.cst, x, ldx, x1, 1);
//             sparse_details_impl::pfspmm_zo(
//                 F, A, blockSize, x1, ldx, y, ldy,
//                 [](Element &a, const Element &b) { a += b; });
//             fflas_delete(x1);
//         }
//     }
// #endif
// }

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    typename Field::Element_ptr y, int ldy,
//                    FieldCategories::ModularTag) {
//     if (A.delayed) {
//         sparse_details::pfspmm(F, A, blockSize, x, ldx, y, ldy,
//                                typename FieldCategories::UnparametricTag(),
//                                typename std::true_type());
//         freduce(F, blockSize, A.m, y, ldy);
//     } else {
//         if (ldx == blockSize && ldy == blockSize) {
//             sparse_details::pfspmm(F, A, blockSize, x, y, A.kmax);
//         } else {
//             sparse_details::pfspmm(F, A, blockSize, x, ldx, y, ldy, A.kmax);
//         }
//     }
// }

} // sparse details

template <class Field, class SM>
inline void fspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                  typename Field::Element_ptr y) {
    sparse_details::init_y(F, A.m, beta, y, typename FieldTraits<Field>::value());
    sparse_details::fspmv<Field, SM>(F, A, x, y, typename FieldTraits<Field>::category(),
                                     typename isZOSparseMatrix<Field, SM>::type());
}

template <class Field, class SM>
inline void pfspmv(const Field &F, const SM &A, typename Field::ConstElement_ptr x, const typename Field::Element &beta,
                   typename Field::Element_ptr y) {
    sparse_details::init_y(F, A.m, beta, y, typename FieldTraits<Field>::value());
    sparse_details::pfspmv<Field, SM>(F, A, x, y, typename FieldTraits<Field>::category(),
                                      typename isZOSparseMatrix<Field, SM>::type());
}

template <class Field, class SM>
inline void fspmm(const Field &F, const SM &A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
                  const typename Field::Element &beta, typename Field::Element_ptr y, int ldy) {
    sparse_details::init_y(F, A.m, blockSize, beta, y, ldy, typename FieldTraits<Field>::value());
    sparse_details::fspmm<Field, SM>(F, A, blockSize, x, ldx, y, ldy, typename FieldTraits<Field>::category(),
                                     typename isZOSparseMatrix<Field, SM>::type());
}

// template <class Field, class SM>
// inline void pfspmm(const Field &F, const SM &A, int blockSize,
//                    typename Field::ConstElement_ptr x, int ldx,
//                    const typename Field::Element &beta,
//                    typename Field::Element_ptr y, int ldy) {
//     sparse_details::init_y(F, A.m, blockSize, beta, y, ldy,
//                            typename FieldTraits<Field>::value());
//     sparse_details::pfspmm<Field, SM>(
//         F, A, blockSize, x, ldx, y, ldy,
//         typename FieldTraits<Field>::category(),
//         typename isZOSparseMatrix<Field, SM>::type());
// }
}

#endif // __FFLASFFPACK_fflas_fflas_sparse_INL
