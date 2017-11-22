/*
 * Copyright (C) 2013 FFLAS-FFPACK group.
 *
 * Extirpé form a m4 macro by Brice Boyer (briceboyer) <boyer.brice@gmail.com>.
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
 *
 */

#define __FFLASFFPACK_CONFIGURATION
#include "fflas-ffpack/config-blas.h"
extern "C" {
	void dgemm_ (const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}
    int main ()
{
	double a[4] = {1.,3.,2.,4.};
	double b[4] = {4.,2.,3.,1.};
	double c[4];
        int n = 2;
        double alpha=1.;
        double beta=0.;
        char tr='N';
	dgemm_(&tr,&tr,&n, &n, &n,&alpha, a,&n,b,&n,&beta,c,&n);
	if ( (c[0]!=8.) && (c[1]!=5.) && (c[2]!=20.) && (c[3]!=13))
		return -1;
	else
		return 0;
}
