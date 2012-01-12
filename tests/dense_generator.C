/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//
/*
 * Coypright (c) FFLAS-FFPACK
 * Written by Clement Pernet <clement.pernet@imag.fr>
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
 */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>



template<class T>
T& myrand (T& r, long size)
{
	  if (size < 0)
		     return r = T( (lrand48() % (-size-size)) + size );
	    else
		       return r = T(  lrand48() % size ) ;
};


int main(int argc, char ** argv)
{

  long ni=10,nj=10,max=100;
  int offset = 0;

	 if (argc > ++offset)
	          ni = atoi( argv[offset] );
	 if (argc > ++offset)
	       nj = atoi( argv[offset] );
	 if (argc > ++offset)
	       max = atoi( argv[offset] );

	 long tmp;
	 printf("%ld %ld M\n", ni, nj);
	 for (long i = 0; i < ni; ++i)
	   for (long j = 0; j < nj; ++j){
	     printf("%ld %ld %ld\n", i+1, j+1,  myrand(tmp, max));
	   }

	 printf("0 0 0\n");

	 return 0;
}
