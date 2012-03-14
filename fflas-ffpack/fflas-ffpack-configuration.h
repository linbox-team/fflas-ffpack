/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2012 FFLAS-FFPACK
 * Written by <brice.boyer@imag.fr>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file fflas-ffpack/fflas-ffpack-configuration.h
 * @ingroup optimise
 * @brief Defaults for optimised values.
 * While \c fflas-ffpack-optimise.h is created by \c configure script,
 * (either left blank or filled by optimiser), this file produces the
 * defaults for the optimised values. If \c fflas-ffpack-optimise.h is not
 * empty, then its values preceeds the defaults here.
 */


#ifndef __FFLASFFPACK_fflas_ffpack_configuration_H
#define __FFLASFFPACK_fflas_ffpack_configuration_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas-ffpack-optimise.h"

// winograd algorithm threshold (for double)
#ifndef __FFLASFFPACK_WINOTHRESHOLD
#define __FFLASFFPACK_WINOTHRESHOLD 1000
#endif

#endif // __FFLASFFPACK_fflas_ffpack_configuration_H
