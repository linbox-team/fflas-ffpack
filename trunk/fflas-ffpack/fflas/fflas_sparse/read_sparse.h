/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla <bastien.vialla@lirmm.fr>
 *              BB <brice.boyer@polsys.lip6.fr>
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

/** @file fflas/fflas_sparse/read_sparse.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_read_sparse_H
#define __FFLASFFPACK_fflas_fflas_sparse_read_sparse_H

// #include <sstream>
// #include <iostream>
#include <fstream> /*  getline */
// #include <string>
// #include <cstdio>
// #include <cstdlib>
#include <iterator> /*  istream_iterator */



namespace FFLAS { namespace details_spmv {
	template <class Field> struct Coo {
	private:
		using Self = Coo<Field>;

	public:
		typename Field::Element val = 0;
		index_t col = 0;
		index_t row = 0;

		Coo() = default;
		Coo(typename Field::Element v, index_t r, index_t c) : val(v), col(c), row(r) {}
		Coo(const Self &) = default;
		Coo(Self &&) = default;

		Self &operator=(const Self &) = default;
		Self &operator=(Self &&) = default;
	};
} // details_spmv
} // FFLAS

namespace FFLAS {

	template <class Field, bool sorted=true, bool read_integer = false>
	void readSmsFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
			   typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz)
	{
		using namespace details_spmv;
		std::ifstream file(path, std::ios::in);
		std::string line;
		std::vector<Coo<Field>> data;
		while (true) { /* comments ? */
			std::getline(file,line);
			if (line.empty()) {
				continue ;
			}
			std::istringstream ligne (line);
			std::string comm ;
			if (ligne >> comm ) {
				if (comm[0] != '%') {
					break;
				}
			}
			else {
				std::cerr << " the impossible happened, continuing for now " << std::endl;
				break;
			}
		}
		bool sms = false ;
		std::istringstream ligne (line);
		std::string nbnz_s ;
		if (ligne >> rowdim >> coldim >> nbnz_s) {
			if (nbnz_s == "M") {
				sms = true ;
				// nnz = 0;
			}
			else
				nnz = std::strtoull(nbnz_s.c_str(),NULL,0);
		}
		else {
			std::cerr << "file " << path << " is not in sms/smf format " << line << std::endl;
			exit(1);
		}

		row = fflas_new<index_t>(rowdim+1);
		std::memset(row,0, sizeof(index_t)*(rowdim+1));
		assert(!row[0] && !row[rowdim]);
		std::vector<index_t> colid((sms)?0:nnz);
		std::vector<typename Field::Element> dat((sms)?0:nnz);


		/* got header */
		if (!rowdim || !coldim)
			exit(-1) ;
		if (!sms && !nnz)
			exit(-1) ;

		size_t i=0,l,c ;
		int64_t d ;
		while (true) {
			if ((!sms) && (i == nnz)){
				break ;
			}
			std::getline(file,line);
			// std::cout << i << ',' << nnz << std::endl;
			if (file.bad() || file.eof())
				exit(-3);
			if (line.empty()){
				continue;
			}

			std::istringstream lign (line);
			if (lign >> l >> c >> d){
				// std::cout << l << ' ' << c << ' ' << d << std::endl;
				if (sms) {
					if (l == 0 && c == 0 && d == 0)
						break ;
					// nnz ++;
				}
				typename Field::Element v;

				assert(l && c);
				f.init(v, d);
				if (!f.isZero(v)) {
					if (!sorted) {
						data.emplace_back(v, l-1, c-1);
					}
					else {
						row[l] += 1 ;
						if (!sms) {
							colid[i] = c-1 ;
							dat[i] = v ;
						}
						else {
							colid.push_back(c-1);
							dat.push_back(v);
						}
					}
				}
			}
			else {
				exit(1);
			}
			++i ;
		}
		if (sms) {
			nnz=dat.size();
		}
		assert(i == nnz);

		col = fflas_new<index_t>(nnz);
		val = fflas_new(f, nnz, 1);

		if (!sorted) {
			assert(nnz == data.size());

			std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
				  return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
				  });
			auto rowmax = (std::max_element(data.begin(), data.end(),
							[](const Coo<Field> &a, const Coo<Field> &b) { return a.row < b.row; }))->row;
			if (rowdim != rowmax + 1) {
				std::cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax << std::endl;
				rowdim = rowmax;
			}


			for (size_t i = 0, end = data.size(); i < end; ++i) {
				val[i] = data[i].val;
				col[i] = data[i].col;
				row[data[i].row+1]+=1;
			}
		}
		else {
			assert(nnz==dat.size());
			for (size_t i = 0, end = nnz; i < end; ++i) {
				val[i] = dat[i];
				col[i] = colid[i];
			}

		}

		for (size_t i = 0, end = rowdim ; i < end; ++i) {
			row[i+1] += row[i] ;
		}
	}

	template <class Field>
	void readSprFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
			   typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz) {
		using namespace details_spmv;
		std::ifstream file(path, std::ios::in);
		std::vector<std::string> tokens;
		std::string line;
		// while(std::getline(file, line) && line.size()!=0);
		std::getline(file, line);
		std::istringstream is(line);
		// std::cout << "line : " << line << std::endl;
		std::copy(std::istream_iterator<std::string>(is), std::istream_iterator<std::string>(),
			  std::back_inserter<std::vector<std::string>>(tokens));
		// std::cout << tokens.size() << std::endl;
		// std::cout << " " << std::stoull(tokens[0]) << " " << std::stoull(tokens[1]) << std::endl;
		rowdim = static_cast<index_t>(std::stoull(tokens[0]));
		coldim = static_cast<index_t>(std::stoull(tokens[1]));
		std::vector<Coo<Field>> data;
		nnz = 0;
		uint64_t itLine = 0;
		while (std::getline(file, line)) {
			tokens.resize(0);
			std::istringstream iss(line);

			std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				  std::back_inserter<std::vector<std::string>>(tokens));

			// if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
			uint64_t nElements = stoull(tokens[0]);
			for (uint64_t i = 0; i < nElements; ++i) {
				index_t c = std::stoull(tokens[2 * i + 1]) - 1;
				typename Field::Element v;
				int64_t vtmp = std::stoll(tokens[2 * (i + 1)]);
				f.init(v, vtmp);
				data.emplace_back(v, itLine, c);
			}
			// typename Field::Element v;
			// f.init(v, std::stol(tokens[2]));
			// index_t r = (index_t)(std::stoull(tokens[0])) - 1;
			// index_t c = (index_t)(std::stoull(tokens[1])) - 1;
			// data.emplace_back(v, r, c);
			// }
			++itLine;
		}
		std::sort(data.begin(), data.end(), [](const Coo<Field> &a, const Coo<Field> &b) {
			  return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
			  });
		auto rowmax = (std::max_element(data.begin(), data.end(),
						[](const Coo<Field> &a, const Coo<Field> &b) { return a.row < b.row; }))->row;
		if (rowdim != rowmax + 1) {
			std::cout << "Matrix row dimension change : " << rowdim << " -> " << rowmax << std::endl;
			rowdim = rowmax;
		}
		row = fflas_new<index_t>(data.size());
		col = fflas_new<index_t>(data.size());
		val = fflas_new(f, data.size(), 1);
		nnz = data.size();
		std::cout << "nnz : " << nnz << std::endl;
		for (size_t i = 0, end = data.size(); i < end; ++i) {
			val[i] = data[i].val;
			col[i] = data[i].col;
			row[i] = data[i].row;
		}
	}

}// FFLAS

#endif /*  __FFLASFFPACK_fflas_fflas_sparse_read_sparse_H */

