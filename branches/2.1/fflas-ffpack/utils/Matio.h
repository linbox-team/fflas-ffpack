/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox,FFLAS-FFPACK
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

#ifndef __FFLASFFPACK_matio_H
#define __FFLASFFPACK_matio_H

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
//#include "fflas-ffpack/fflas/fflas.h"
#include "fflas_memory.h"

// Reading and writing matrices over field

// Reading a matrice from a (eventually zipped) file
template<class Field>
typename Field::Element_ptr read_field(const Field& F, const char * mat_file,int* tni,int* tnj)
{
	char *UT = NULL;
	const char* File_Name;
	int is_gzipped = 0;
	size_t s = strlen(mat_file);
	typename Field::Element_ptr X = NULL;
	if ((mat_file[--s] == 'z') &&
	    (mat_file[--s] == 'g') &&
	    (mat_file[--s] == '.')) {
		is_gzipped = 1;
		char tmp_nam[] = "/tmp/bbXXXXXX_";
		if (mkstemp(tmp_nam))
			printf("Error opening file]\n");
		File_Name  = tmp_nam;

		UT = FFLAS::fflas_new<char>(s+34+strlen(File_Name));
		sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
		if (system(UT))
			printf("Error uncompressing file\n");
		sprintf(UT,"\\rm %s", File_Name);
	} else {
		File_Name = mat_file;
	}

	FILE* FileDes = fopen(File_Name, "r");
	if (FileDes != NULL) {
		char  tmp [200];// unsigned long tni, tnj;
		if (fscanf(FileDes,"%d %d %199s\n",tni, tnj, tmp)<0)
			printf("Error Reading first line of file \n");
		int n=*tni;
		int p=*tnj;
		X = FFLAS::fflas_new<typename Field::Element>(n*p);
		for (int i=0;i<n*p;++i)
			F.assign(X[i], F.zero);
		long i,j; long val;
		if(fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val)<0)
			printf("Read Error\n");
		while(i && j) {
			F.init(X[p*(i-1)+j-1],val);
			if(fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val)<0)
				printf("Read Error\n");
		}
		fclose(FileDes);
	}

	if (is_gzipped)
		if (system(UT))
			printf("Error uncompressing file\n");
	if (UT != NULL)
		FFLAS::fflas_delete( UT);
	return X;
}

// Displays a matrix
template<class Field>
std::ostream& write_field(const Field& F,std::ostream& c,
			  typename Field::ConstElement_ptr E,
			  int n, int m, int id, bool mapleFormat = false, bool column_major=false)
{

	    //typename Field::Element tmp;
	// double tmp;
//	Givaro::Integer tmp;
	typename Field::Element tmp;
	F.init(tmp);
	if (mapleFormat) c << "Matrix(" << n <<',' << m << ",\n[" ;
	for (int i = 0; i<n;++i){
		if (mapleFormat) c << '[';
		for (int j=0; j<m;++j){
			if (column_major)
				    //F.convert(tmp,*(E+i+id*j));
				    tmp = *(E+i+id*j);
				
			else
//				F.convert(tmp,*(E+j+id*i));
				tmp =*(E+j+id*i);
			c << tmp;
			if (mapleFormat && j<m-1) c << ',';
			c << ' ';
		}
		if (mapleFormat) c << ']';
		if (mapleFormat && i<n-1) c << ',';
		if (!mapleFormat) c << std::endl;
	}
	if (mapleFormat) c << "]);";
	return c ;
}

inline std::ostream& write_perm (std::ostream& c, const size_t* P, size_t N){
	c<<"[ ";
	for (size_t i=0; i<N; ++i)
		c<<P[i]<<" ";
	c<<"]"<<std::endl;
	return c;
}
#endif //__FFLASFFPACK_matio_H
