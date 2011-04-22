/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* field/modular-balanced-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005 Clement Pernet
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * and Clement Pernet <Clement.Pernet@imag.fr>
 * and Brice Boyer <bboyer@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/*! @file field/modular-balanced-double.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c double .
 */

#ifndef __FFLAFLAS_modular_balanced_H
#define __FFLAFLAS_modular_balanced_H

#include <math.h>
#include "fflas-ffpack/field/modular-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"

namespace FFPACK {

template <class Element>
class ModularBalanced;

}


#include "fflas-ffpack/field/modular-balanced-double.h"
#include "fflas-ffpack/field/modular-balanced-float.h"
#include "fflas-ffpack/field/modular-balanced-int32.h"
#include "fflas-ffpack/field/modular-balanced-int64.h"

#endif // __FFLAFLAS_modular_balanced_double_H

