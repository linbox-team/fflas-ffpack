/* fflas_L1_inst.h
 * Copyright (C) 2015 FFLAS-FFPACK group
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */
#ifndef __FFLAS_L1_INST_H
#define __FFLAS_L1_INST_H

#include "givaro/modular.h"
#include "givaro/modular-balanced.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_helpers.inl"

#ifdef INST_OR_DECL
#undef INST_OR_DECL
#endif
#define INST_OR_DECL <>

#define FFLAS_FIELD Givaro::ModularBalanced
#define FFLAS_ELT double
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#define FFLAS_ELT float
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#define FFLAS_ELT int64_t
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#undef FFLAS_FIELD

#define FFLAS_FIELD Givaro::Modular
#define FFLAS_ELT double
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#define FFLAS_ELT float
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#define FFLAS_ELT int64_t
#include "fflas_L1_inst_implem.inl"
#undef FFLAS_ELT
#undef FFLAS_FIELD

#endif //__FFLAS_L1_INST_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
