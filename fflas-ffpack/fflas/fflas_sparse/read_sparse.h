/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla <bastien.vialla@lirmm.fr>
 *              Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#include "fflas-ffpack/fflas-ffpack-config.h"

// #include <sstream>
// #include <iostream>
#include <fstream> /*  getline */
#include <string>
// #include <cstdio>
#include <cstdlib>
#include <iterator> /*  istream_iterator */



namespace FFLAS { namespace details_spmv {
    template <class Field> struct Coo {
    private:
        using Self = Coo<Field>;

    public:
        typename Field::Element val = 0;
        index_t col = 0;
        index_t row = 0;
        bool deleted = false;

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
        val = fflas_new(f, nnz);

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


            for (size_t j = 0, end = data.size(); j < end; ++j) {
                val[j] = data[j].val;
                col[j] = data[j].col;
                row[data[j].row+1]+=1;
            }
        }
        else {
            assert(nnz==dat.size());
            for (size_t j = 0, end = nnz; j < end; ++j) {
                val[j] = dat[j];
                col[j] = colid[j];
            }

        }

        for (size_t j = 0, end = rowdim ; j < end; ++j) {
            row[j+1] += row[j] ;
        }
    }

    template <class Field>
    void readSprFormat(const std::string &path, const Field &f, index_t *&row, index_t *&col,
                       typename Field::Element_ptr &val, index_t &rowdim, index_t &coldim, uint64_t &nnz)
    {
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
        while (!file.eof()) {
            std::getline(file, line);
            tokens.resize(0);
            std::istringstream iss(line);

            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                      std::back_inserter<std::vector<std::string>>(tokens));

            // if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
            uint64_t nElements = std::stoull(tokens[0]);
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
        val = fflas_new(f, data.size());
        nnz = data.size();
        std::cout << "nnz : " << nnz << std::endl;
        for (size_t i = 0, end = data.size(); i < end; ++i) {
            val[i] = data[i].val;
            col[i] = data[i].col;
            row[i] = data[i].row;
        }
    }

#define DNS_BIN_VER 0
#define mask_t uint64_t

    template<class Field, class T>
    struct readMyMachineType {
        typedef typename Field::Element     Element ;
        typedef typename Field::Element_ptr Element_ptr ;
        void operator() (const Field &F,
                         Element & modulo,
                         Element_ptr val,
                         std::ifstream & file,
                         const uint64_t dims,
                         const mask_t data_type,
                         const mask_t field_desc);
    };

    template<class Field>
    struct readMyMachineType<Field,mpz_t> {
        typedef typename Field::Element     Element ;
        typedef typename Field::Element_ptr Element_ptr ;
        void operator() (const Field &F,
                         Element & modulo,
                         Element_ptr val,
                         std::ifstream & file,
                         const uint64_t dims,
                         const mask_t data_type,
                         const mask_t field_desc);
    };



    template<class Field, typename T>
    void readMyMachineType<Field,T>:: operator() (const Field &F,
                                                  Element & modulo,
                                                  Element_ptr val,
                                                  std::ifstream & file,
                                                  const uint64_t dims,
                                                  const mask_t data_type,
                                                  const mask_t field_desc)
    {
        if (field_desc ==1) { /*  modulo */
            T modulo_read ;
            file.read((char*) &modulo_read, sizeof(T));
            F.init(modulo, modulo_read);
        }
        /*  do something with field_desc and multiprec... */
        T * data_read = fflas_new<T>(dims);
        file.read((char*)data_read,sizeof(T));
        /* TODO freduce ? */
        for (size_t i = 0 ; i< dims ; ++i) {
            F.init(val[i], data_read[i]);
        }
    }

    template<class Field>
    void readMyMachineType<Field,mpz_t>:: operator() (const Field &F,
                                                      typename Field::Element & modulo,
                                                      typename Field::Element_ptr val,
                                                      std::ifstream & file,
                                                      const uint64_t dims,
                                                      const mask_t data_type,
                                                      const mask_t field_desc)
    {
        /* need to use FILE * instead of std::ifstream */
        throw("not implemented, use mpz_in_raw, but FILE*...");
    }


    template<class T>
    std::enable_if<std::is_integral<T>::value,int>
    getDataType()
    {
        return (1<<(sizeof(T)-1))+ std::is_unsigned<T>::value ;
    }

    template<class T>
    std::enable_if<std::is_floating_point<T>::value,int>
    getDataType()
    {
        return (1<<8)+std::is_same<T,double>::value ;
    }

    template<class T>
    std::enable_if<std::is_same<T,mpz_t>::value,int>
    getDataType()
    {
        return (1<<16) ;
    }

    template<class T>
    int getDataType()
    {
        return -1 ;
    }


    template<class Field>
    void readMachineType(const Field &F,
                         typename Field::Element & modulo,
                         typename Field::Element_ptr val,
                         std::ifstream & file,
                         const uint64_t dims,
                         const mask_t data_type,
                         const mask_t field_desc)
    {
        // switch(data_type) {
        // case (1<<0) + 0 :
        // 	readMyMachineType<Field,int8_t   >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<0) + 1 :
        // 	readMyMachineType<Field,uint8_t  >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<1) + 0 :
        // 	readMyMachineType<Field,int16_t  >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<1) + 1 :
        // 	readMyMachineType<Field,uint16_t >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<2) + 0 :
        // 	readMyMachineType<Field,int32_t  >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<2) + 0 :
        // 	readMyMachineType<Field,uint32_t >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<3) + 0 :
        // 	readMyMachineType<Field,int64_t  >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<3) + 0 :
        // 	readMyMachineType<Field,uint64_t >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<8) :
        // 	readMyMachineType<Field,float    >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<8)+1 :
        // 	readMyMachineType<Field,double   >() (F,val, modulo, file,dims,data_type,field_desc);
        // case (1<<16) :
        // 	readMyMachineType<Field,mpz_t    >() (F,val, modulo, file,dims,data_type,field_desc);
        // default :
        // 	throw("bad data type descriptor");
        // }

    }

    template<class Field>
    void readDnsFormat(const std::string &path, const Field &F, index_t &rowdim, index_t &coldim,
                       typename Field::Element_ptr &val)
    {
        std::ifstream file(path, std::ifstream::binary);
        mask_t magic, field_desc, data_type ;
        typename Field::Element modulo ;

        file.read((char*) &magic     , sizeof(int64_t)) ;
        if (magic != DNS_BIN_VER) {
            throw("bad version");
        }
        file.read((char*) &field_desc, sizeof(int64_t)) ;
        file.read((char*) &data_type , sizeof(int64_t)) ;
        file.read((char*) &rowdim , sizeof(int64_t)) ;
        file.read((char*) &coldim , sizeof(int64_t)) ;
        val = fflas_new(F,rowdim,coldim);
        readMachineType(F,val, modulo, file,rowdim*coldim,field_desc,data_type);

    }

    template<class Field>
    void writeDnsFormat(const std::string &path, const Field &F, const index_t &rowdim, const index_t &coldim,
                        typename Field::Element_ptr A, index_t ldA)
    {
        typedef typename Field::Element Element ;
        std::ofstream file(path, std::ofstream::binary);
        mask_t field_desc = getFieldDesc(F);
        mask_t magic = DNS_BIN_VER ;
        mask_t data_type = getDataType<Element>(F);
        Element modulo ;

        file.write((char*) &magic     , sizeof(int64_t)) ;
        file.write((char*) &field_desc, sizeof(int64_t)) ;
        file.write((char*) &data_type , sizeof(int64_t)) ;
        file.write((char*) &rowdim , sizeof(int64_t)) ;
        file.write((char*) &coldim , sizeof(int64_t)) ;
        // writeMachineType(F,A, modulo, file,rowdim,coldim,lda,field_desc,data_type);

    }

}// FFLAS

#endif /*  __FFLASFFPACK_fflas_fflas_sparse_read_sparse_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
