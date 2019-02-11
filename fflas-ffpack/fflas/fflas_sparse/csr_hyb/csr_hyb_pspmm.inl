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

#ifndef __FFLASFFPACK_fflas_sparse_CSR_HYB_pspmm_INL
#define __FFLASFFPACK_fflas_sparse_CSR_HYB_pspmm_INL

namespace FFLAS {
    namespace sparse_details_impl {
        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, typename Field::Element_ptr y, FieldCategories::GenericTag) {
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              F.subin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                              F.subin(y[i * blockSize + k + 1], x[A.col[j] * blockSize + k + 1]);
                              F.subin(y[i * blockSize + k + 2], x[A.col[j] * blockSize + k + 2]);
                              F.subin(y[i * blockSize + k + 3], x[A.col[j] * blockSize + k + 3]);
                              }
                              for (; k < blockSize; ++k)
                              F.subin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              F.addin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                              F.addin(y[i * blockSize + k + 1], x[A.col[j] * blockSize + k + 1]);
                              F.addin(y[i * blockSize + k + 2], x[A.col[j] * blockSize + k + 2]);
                              F.addin(y[i * blockSize + k + 3], x[A.col[j] * blockSize + k + 3]);
                              }
                              for (; k < blockSize; ++k)
                                  F.addin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                                      F.axpyin(y[i * blockSize + k], A.dat[startDat + k], x[A.col[j] * blockSize + k]);
                                      F.axpyin(y[i * blockSize + k + 1], A.dat[startDat + k], x[A.col[j] * blockSize + k + 1]);
                                      F.axpyin(y[i * blockSize + k + 2], A.dat[startDat + k], x[A.col[j] * blockSize + k + 2]);
                                      F.axpyin(y[i * blockSize + k + 3], A.dat[startDat + k], x[A.col[j] * blockSize + k + 3]);
                                  }
                                  for (; k < blockSize; ++k)
                                      F.axpyin(y[i * blockSize + k], A.dat[startDat + k], x[A.col[j] * blockSize + k]);
                              }
                              }
                              });
#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.subin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                        F.subin(y[i * blockSize + k + 1], x[A.col[j] * blockSize + k + 1]);
                        F.subin(y[i * blockSize + k + 2], x[A.col[j] * blockSize + k + 2]);
                        F.subin(y[i * blockSize + k + 3], x[A.col[j] * blockSize + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.subin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.addin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                        F.addin(y[i * blockSize + k + 1], x[A.col[j] * blockSize + k + 1]);
                        F.addin(y[i * blockSize + k + 2], x[A.col[j] * blockSize + k + 2]);
                        F.addin(y[i * blockSize + k + 3], x[A.col[j] * blockSize + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.addin(y[i * blockSize + k], x[A.col[j] * blockSize + k]);
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.axpyin(y[i * blockSize + k], A.dat[startDat + k], x[A.col[j] * blockSize + k]);
                        F.axpyin(y[i * blockSize + k + 1], A.dat[startDat + k], x[A.col[j] * blockSize + k + 1]);
                        F.axpyin(y[i * blockSize + k + 2], A.dat[startDat + k], x[A.col[j] * blockSize + k + 2]);
                        F.axpyin(y[i * blockSize + k + 3], A.dat[startDat + k], x[A.col[j] * blockSize + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.axpyin(y[i * blockSize + k], A.dat[startDat + k], x[A.col[j] * blockSize + k]);
                }
            }
#endif
        }

        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, int ldx, typename Field::Element_ptr y, int ldy,
                           FieldCategories::GenericTag) {
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize, ldx, ldy](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              F.subin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                              F.subin(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                              F.subin(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                              F.subin(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
                              }
                              for (; k < blockSize; ++k)
                              F.subin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              F.addin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                              F.addin(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                              F.addin(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                              F.addin(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
                              }
                              for (; k < blockSize; ++k)
                                  F.addin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                                      F.axpyin(y[i * ldy + k], A.dat[startDat + k], x[A.col[j] * ldx + k]);
                                      F.axpyin(y[i * ldy + k + 1], A.dat[startDat + k], x[A.col[j] * ldx + k + 1]);
                                      F.axpyin(y[i * ldy + k + 2], A.dat[startDat + k], x[A.col[j] * ldx + k + 2]);
                                      F.axpyin(y[i * ldy + k + 3], A.dat[startDat + k], x[A.col[j] * ldx + k + 3]);
                                  }
                                  for (; k < blockSize; ++k)
                                      F.axpyin(y[i * ldy + k], A.dat[startDat + k], x[A.col[j] * ldx + k]);
                              }
                              }
                              });
#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.subin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                        F.subin(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                        F.subin(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                        F.subin(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.subin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.addin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                        F.addin(y[i * ldy + k + 1], x[A.col[j] * ldx + k + 1]);
                        F.addin(y[i * ldy + k + 2], x[A.col[j] * ldx + k + 2]);
                        F.addin(y[i * ldy + k + 3], x[A.col[j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.addin(y[i * ldy + k], x[A.col[j] * ldx + k]);
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        F.axpyin(y[i * ldy + k], A.dat[startDat + k], x[A.col[j] * ldx + k]);
                        F.axpyin(y[i * ldy + k + 1], A.dat[startDat + k], x[A.col[j] * ldx + k + 1]);
                        F.axpyin(y[i * ldy + k + 2], A.dat[startDat + k], x[A.col[j] * ldx + k + 2]);
                        F.axpyin(y[i * ldy + k + 3], A.dat[startDat + k], x[A.col[j] * ldx + k + 3]);
                    }
                    for (; k < blockSize; ++k)
                        F.axpyin(y[i * ldy + k], A.dat[startDat + k], x[A.col[j] * ldx + k]);
                }
            }
#endif
        }

        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, typename Field::Element_ptr y,
                           FieldCategories::UnparametricTag) {
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                              y[i * blockSize + k + 1] -= x[A.col[j] * blockSize + k + 1];
                              y[i * blockSize + k + 2] -= x[A.col[j] * blockSize + k + 2];
                              y[i * blockSize + k + 3] -= x[A.col[j] * blockSize + k + 3];
                              }
                              for (; k < blockSize; ++k)
                              y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                              y[i * blockSize + k + 1] += x[A.col[j] * blockSize + k + 1];
                              y[i * blockSize + k + 2] += x[A.col[j] * blockSize + k + 2];
                              y[i * blockSize + k + 3] += x[A.col[j] * blockSize + k + 3];
                              }
                              for (; k < blockSize; ++k)
                                  y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                                      y[i * blockSize + k] += A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                                      y[i * blockSize + k + 1] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 1];
                                      y[i * blockSize + k + 2] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 2];
                                      y[i * blockSize + k + 3] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 3];
                                  }
                                  for (; k < blockSize; ++k)
                                      y[i * blockSize + k] += A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                              }
                              }
                              });
#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                        y[i * blockSize + k + 1] -= x[A.col[j] * blockSize + k + 1];
                        y[i * blockSize + k + 2] -= x[A.col[j] * blockSize + k + 2];
                        y[i * blockSize + k + 3] -= x[A.col[j] * blockSize + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                        y[i * blockSize + k + 1] += x[A.col[j] * blockSize + k + 1];
                        y[i * blockSize + k + 2] += x[A.col[j] * blockSize + k + 2];
                        y[i * blockSize + k + 3] += x[A.col[j] * blockSize + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * blockSize + k] += A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                        y[i * blockSize + k + 1] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 1];
                        y[i * blockSize + k + 2] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 2];
                        y[i * blockSize + k + 3] += A.dat[startDat + j] * x[A.col[j] * blockSize + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] += A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                }
            }
#endif
        }

        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, int ldx, typename Field::Element_ptr y, int ldy,
                           FieldCategories::UnparametricTag) {
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize, ldx, ldy](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              y[i * ldy + k] -= x[A.col[j] * ldx + k];
                              y[i * ldy + k + 1] -= x[A.col[j] * ldx + k + 1];
                              y[i * ldy + k + 2] -= x[A.col[j] * ldx + k + 2];
                              y[i * ldy + k + 3] -= x[A.col[j] * ldx + k + 3];
                              }
                              for (; k < blockSize; ++k)
                              y[i * ldy + k] -= x[A.col[j] * ldx + k];
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                              y[i * ldy + k] += x[A.col[j] * ldx + k];
                              y[i * ldy + k + 1] += x[A.col[j] * ldx + k + 1];
                              y[i * ldy + k + 2] += x[A.col[j] * ldx + k + 2];
                              y[i * ldy + k + 3] += x[A.col[j] * ldx + k + 3];
                              }
                              for (; k < blockSize; ++k)
                                  y[i * ldy + k] += x[A.col[j] * ldx + k];
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                                      y[i * ldy + k] += A.dat[startDat + j] * x[A.col[j] * ldx + k];
                                      y[i * ldy + k + 1] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 1];
                                      y[i * ldy + k + 2] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 2];
                                      y[i * ldy + k + 3] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 3];
                                  }
                                  for (; k < blockSize; ++k)
                                      y[i * ldy + k] += A.dat[startDat + j] * x[A.col[j] * ldx + k];
                              }
                              }
                              });
#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] -= x[A.col[j] * ldx + k];
                        y[i * ldy + k + 1] -= x[A.col[j] * ldx + k + 1];
                        y[i * ldy + k + 2] -= x[A.col[j] * ldx + k + 2];
                        y[i * ldy + k + 3] -= x[A.col[j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] -= x[A.col[j] * ldx + k];
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] += x[A.col[j] * ldx + k];
                        y[i * ldy + k + 1] += x[A.col[j] * ldx + k + 1];
                        y[i * ldy + k + 2] += x[A.col[j] * ldx + k + 2];
                        y[i * ldy + k + 3] += x[A.col[j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += x[A.col[j] * ldx + k];
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 4); k += 4) {
                        y[i * ldy + k] += A.dat[startDat + j] * x[A.col[j] * ldx + k];
                        y[i * ldy + k + 1] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 1];
                        y[i * ldy + k + 2] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 2];
                        y[i * ldy + k + 3] += A.dat[startDat + j] * x[A.col[j] * ldx + k + 3];
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += A.dat[startDat + j] * x[A.col[j] * ldx + k];
                }
            }
#endif
        }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        template <class Field, class LFunc, class SFunc>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, typename Field::Element_ptr y, LFunc &&lfunc, SFunc &&sfunc,
                           FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize, lfunc, sfunc](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              vect_t vx1, vx2, vy1, vy2, vdat;
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                              vy1 = lfunc(y + i * blockSize + k);
                              vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                              vx1 = lfunc(y + A.col[j] * blockSize + k);
                              vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                              sfunc(y + i * blockSize + k, simd::sub(vy1, vx1));
                              sfunc(y + i * blockSize + k + simd::vect_size, simd::sub(vy2, vx2));
                              }
                              for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                              vy1 = lfunc(y + i * blockSize + k);
                              vx1 = lfunc(y + A.col[j] * blockSize + k);
                              sfunc(y + i * blockSize + k, simd::sub(vy1, vx1));
                              }
                              for (; k < blockSize; ++k)
                                  y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                                      vy1 = lfunc(y + i * blockSize + k);
                                      vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                                      vx1 = lfunc(y + A.col[j] * blockSize + k);
                                      vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                                      sfunc(y + i * blockSize + k, simd::add(vy1, vx1));
                                      sfunc(y + i * blockSize + k + simd::vect_size, simd::add(vy2, vx2));
                                  }
                                  for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                                      vy1 = lfunc(y + i * blockSize + k);
                                      vx1 = lfunc(y + A.col[j] * blockSize + k);
                                      sfunc(y + i * blockSize + k, simd::add(vy1, vx1));
                                  }
                                  for (; k < blockSize; ++k)
                                      y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  for (uint64_t j = start; j < stop; ++j) {
                                      size_t k = 0;
                                      vdat = simd::set1(A.dat[startDat + j]);
                                      for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                                          vy1 = lfunc(y + i * blockSize + k);
                                          vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                                          vx1 = lfunc(y + A.col[j] * blockSize + k);
                                          vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                                          sfunc(y + i * blockSize + k, simd::fmadd(vy1, vdat, vx1));
                                          sfunc(y + i * blockSize + k + simd::vect_size, simd::fmadd(vy2, vdat, vx2));
                                      }
                                      for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                                          vy1 = lfunc(y + i * blockSize + k);
                                          vx1 = lfunc(y + A.col[j] * blockSize + k);
                                          sfunc(y + i * blockSize + k, simd::fmadd(vy1, vdat, vx1));
                                      }
                                      for (; k < blockSize; ++k)
                                          y[i * blockSize + k] -= A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                                  }
                              }
                              }
                              });
#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                vect_t vx1, vx2, vy1, vy2, vdat;
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = lfunc(y + i * blockSize + k);
                        vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                        vx1 = lfunc(y + A.col[j] * blockSize + k);
                        vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                        sfunc(y + i * blockSize + k, simd::sub(vy1, vx1));
                        sfunc(y + i * blockSize + k + simd::vect_size, simd::sub(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = lfunc(y + i * blockSize + k);
                        vx1 = lfunc(y + A.col[j] * blockSize + k);
                        sfunc(y + i * blockSize + k, simd::sub(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] -= x[A.col[j] * blockSize + k];
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = lfunc(y + i * blockSize + k);
                        vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                        vx1 = lfunc(y + A.col[j] * blockSize + k);
                        vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                        sfunc(y + i * blockSize + k, simd::add(vy1, vx1));
                        sfunc(y + i * blockSize + k + simd::vect_size, simd::add(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = lfunc(y + i * blockSize + k);
                        vx1 = lfunc(y + A.col[j] * blockSize + k);
                        sfunc(y + i * blockSize + k, simd::add(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * blockSize + k] += x[A.col[j] * blockSize + k];
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    for (uint64_t j = start; j < stop; ++j) {
                        size_t k = 0;
                        vdat = simd::set1(A.dat[startDat + j]);
                        for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                            vy1 = lfunc(y + i * blockSize + k);
                            vy2 = lfunc(y + i * blockSize + k + simd::vect_size);
                            vx1 = lfunc(y + A.col[j] * blockSize + k);
                            vx2 = lfunc(y + A.col[j] * blockSize + k + simd::vect_size);
                            sfunc(y + i * blockSize + k, simd::fmadd(vy1, vdat, vx1));
                            sfunc(y + i * blockSize + k + simd::vect_size, simd::fmadd(vy2, vdat, vx2));
                        }
                        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                            vy1 = lfunc(y + i * blockSize + k);
                            vx1 = lfunc(y + A.col[j] * blockSize + k);
                            sfunc(y + i * blockSize + k, simd::fmadd(vy1, vdat, vx1));
                        }
                        for (; k < blockSize; ++k)
                            y[i * blockSize + k] -= A.dat[startDat + j] * x[A.col[j] * blockSize + k];
                    }
                }
            }
#endif
        }

        template <class Field, class LFunc, class SFunc>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, int ldx, typename Field::Element_ptr y, int ldy, LFunc &&lfunc,
                           SFunc &&sfunc, FieldCategories::UnparametricTag) {
            using simd = Simd<typename Field::Element>;
            using vect_t = typename simd::vect_t;
#ifdef __FFLASFFPACK_USE_TBB
            tbb::parallel_for(tbb::blocked_range<index_t>(0, A.m),
                              [&F, &A, &x, &y, blockSize, ldx, ldy, lfunc, sfunc](const tbb::blocked_range<index_t> &r) {
                              for (index_t i = r.begin(), end = r.end(); i < end; ++i) {
                              vect_t vx1, vx2, vy1, vy2, vdat;
                              index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                              for (uint64_t j = start; j < stop; ++j) {
                              size_t k = 0;
                              for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                              vy1 = lfunc(y + i * ldy + k);
                              vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                              vx1 = lfunc(y + A.col[j] * ldx + k);
                              vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                              sfunc(y + i * ldy + k, simd::sub(vy1, vx1));
                              sfunc(y + i * ldy + k + simd::vect_size, simd::sub(vy2, vx2));
                              }
                              for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                              vy1 = lfunc(y + i * ldy + k);
                              vx1 = lfunc(y + A.col[j] * ldx + k);
                              sfunc(y + i * ldy + k, simd::sub(vy1, vx1));
                              }
                              for (; k < blockSize; ++k)
                                  y[i * ldy + k] -= x[A.col[j] * ldx + k];
                              }
                              start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                              for (uint64_t j = start; j < stop; ++j) {
                                  size_t k = 0;
                                  for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                                      vy1 = lfunc(y + i * ldy + k);
                                      vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                                      vx1 = lfunc(y + A.col[j] * ldx + k);
                                      vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                                      sfunc(y + i * ldy + k, simd::add(vy1, vx1));
                                      sfunc(y + i * ldy + k + simd::vect_size, simd::add(vy2, vx2));
                                  }
                                  for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                                      vy1 = lfunc(y + i * ldy + k);
                                      vx1 = lfunc(y + A.col[j] * ldx + k);
                                      sfunc(y + i * ldy + k, simd::add(vy1, vx1));
                                  }
                                  for (; k < blockSize; ++k)
                                      y[i * ldy + k] += x[A.col[j] * ldx + k];
                              }
                              start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                              index_t startDat = A.st[4 * i + 3];
                              for (uint64_t j = start; j < stop; ++j) {
                                  for (uint64_t j = start; j < stop; ++j) {
                                      size_t k = 0;
                                      vdat = simd::set1(A.dat[startDat + j]);
                                      for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                                          vy1 = lfunc(y + i * ldy + k);
                                          vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                                          vx1 = lfunc(y + A.col[j] * ldx + k);
                                          vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                                          sfunc(y + i * ldy + k, simd::fmadd(vy1, vdat, vx1));
                                          sfunc(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vdat, vx2));
                                      }
                                      for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                                          vy1 = lfunc(y + i * ldy + k);
                                          vx1 = lfunc(y + A.col[j] * ldx + k);
                                          sfunc(y + i * ldy + k, simd::fmadd(vy1, vdat, vx1));
                                      }
                                      for (; k < blockSize; ++k)
                                          y[i * ldy + k] -= A.dat[startDat + j] * x[A.col[j] * ldx + k];
                                  }
                              }
                              }
                              });

#else
#pragma omp parallel for
            for (uint64_t i = 0; i < A.m; ++i) {
                vect_t vx1, vx2, vy1, vy2, vdat;
                index_t start = A.st[4 * i], stop = A.st[4 * i + 1];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = lfunc(y + i * ldy + k);
                        vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                        vx1 = lfunc(y + A.col[j] * ldx + k);
                        vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                        sfunc(y + i * ldy + k, simd::sub(vy1, vx1));
                        sfunc(y + i * ldy + k + simd::vect_size, simd::sub(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = lfunc(y + i * ldy + k);
                        vx1 = lfunc(y + A.col[j] * ldx + k);
                        sfunc(y + i * ldy + k, simd::sub(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] -= x[A.col[j] * ldx + k];
                }
                start = A.st[4 * i + 1], stop = A.st[4 * i + 2];
                for (uint64_t j = start; j < stop; ++j) {
                    size_t k = 0;
                    for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                        vy1 = lfunc(y + i * ldy + k);
                        vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                        vx1 = lfunc(y + A.col[j] * ldx + k);
                        vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                        sfunc(y + i * ldy + k, simd::add(vy1, vx1));
                        sfunc(y + i * ldy + k + simd::vect_size, simd::add(vy2, vx2));
                    }
                    for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                        vy1 = lfunc(y + i * ldy + k);
                        vx1 = lfunc(y + A.col[j] * ldx + k);
                        sfunc(y + i * ldy + k, simd::add(vy1, vx1));
                    }
                    for (; k < blockSize; ++k)
                        y[i * ldy + k] += x[A.col[j] * ldx + k];
                }
                start = A.st[4 * i + 2], stop = A.st[4 * (i + 1)];
                index_t startDat = A.st[4 * i + 3];
                for (uint64_t j = start; j < stop; ++j) {
                    for (uint64_t j = start; j < stop; ++j) {
                        size_t k = 0;
                        vdat = simd::set1(A.dat[startDat + j]);
                        for (; k < ROUND_DOWN(blockSize, 2 * simd::vect_size); k += 2 * simd::vect_size) {
                            vy1 = lfunc(y + i * ldy + k);
                            vy2 = lfunc(y + i * ldy + k + simd::vect_size);
                            vx1 = lfunc(y + A.col[j] * ldx + k);
                            vx2 = lfunc(y + A.col[j] * ldx + k + simd::vect_size);
                            sfunc(y + i * ldy + k, simd::fmadd(vy1, vdat, vx1));
                            sfunc(y + i * ldy + k + simd::vect_size, simd::fmadd(vy2, vdat, vx2));
                        }
                        for (; k < ROUND_DOWN(blockSize, simd::vect_size); k += simd::vect_size) {
                            vy1 = lfunc(y + i * ldy + k);
                            vx1 = lfunc(y + A.col[j] * ldx + k);
                            sfunc(y + i * ldy + k, simd::fmadd(vy1, vdat, vx1));
                        }
                        for (; k < blockSize; ++k)
                            y[i * ldy + k] -= A.dat[startDat + j] * x[A.col[j] * ldx + k];
                    }
                }
            }
#endif
        }
#endif

        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, typename Field::Element_ptr y, const int64_t kmax) {
            // TODO
        }

        template <class Field>
        inline void pfspmm(const Field &F, const Sparse<Field, SparseMatrix_t::CSR_HYB> &A, size_t blockSize,
                           typename Field::ConstElement_ptr x, int ldx, typename Field::Element_ptr y, int ldy,
                           const int64_t kmax) {
            // TODO
        }

    } // csr_hyb_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_CSR_HYB_pspmm_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
