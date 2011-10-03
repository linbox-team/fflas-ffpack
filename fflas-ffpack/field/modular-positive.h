/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/modular-positive.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FFLASFFPACK_modular_positive_H
#define __FFLASFFPACK_modular_positive_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"

namespace FFPACK {

	template <class Element>
	class Modular;

} // FFPACK

#include "fflas-ffpack/field/modular-float.h"
#include "fflas-ffpack/field/modular-double.h"
#include "fflas-ffpack/field/modular-int32.h"
#ifdef __x86_64__
#include "fflas-ffpack/field/modular-int64.h"
#endif

#endif // __FFLASFFPACK_modular_positive_H
