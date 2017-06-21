/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/print-utils.h
 * Copyright (C) 2011, Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * 					   Bastien Vialla <bastien.vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_print_utils_H
#define __FFLASFFPACK_print_utils_H

#include <fflas-ffpack/fflas-ffpack-config.h>
#include <math.h>
#include <givaro/givprint.h>
#include <givaro/givinteger.h>

namespace std
{


	template<typename T>
	std::ostream& write_matrix(std::ostream& out, Givaro::Integer p, size_t m, size_t n, T* C, size_t ldc){
		
		size_t www(size_t((double(p.bitsize())*log(2.))/log(10.)));
		out<<"Matrix("<<m<<','<<n<<",[[";
		out.width(www+1);
		out<<std::right<<C[0];
		for (size_t j=1;j<n;++j){
			out<<',';
			out.width(www);
			out<<std::right<<C[j];
		}
		out<<']';
		for (size_t i=1;i<m;++i){ 
			out<<",[";
			out.width(www+1);
			out<<std::right<<C[i*ldc];
			for (size_t j=1;j<n;++j){
				out<<',';
				out.width(www);
				out<<std::right<<C[i*ldc+j];
			}
			out<<']';
		}
		return out<<"]);"<<std::endl;
	}
	
	template<typename Field>
	std::ostream& write_matrix(std::ostream& out, const Field& F, size_t m, size_t n, typename Field::ConstElement_ptr C, size_t ldc){
		
		out<<"Matrix("<<m<<','<<n<<",[[";
		F.write(out,C[0]);
		for (size_t j=1;j<n;++j){
			out<<',';
			F.write(out,C[j]);
		}
		out<<']';
		for (size_t i=1;i<m;++i){ 
			out<<",[";
			F.write(out,C[i*ldc]);
			for (size_t j=1;j<n;++j){
				out<<',';
				F.write(out,C[i*ldc+j]);
			}
			out<<']';
		}
		return out<<"]);"<<std::endl;
	}
	
}

#endif // __FFLASFFPACK_print_utils_H
