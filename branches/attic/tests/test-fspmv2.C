/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :
 *        Bastien Vialla <bastien.vialla@lirmm.fr>
 * This file is Free Software and part of FFLAS-FFPACK.
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

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_fspmv.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/field/modular-double.h"
#include "fflas-ffpack/field/unparametric.h"
// #include "fflas-ffpack/utils/flimits.h"
// #include "fllas-ffpack/utils/print-utils.h"

#include "fflas-ffpack/field/rns.h" 

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iterator>

using namespace FFLAS;
using namespace FFPACK;
using namespace std;

namespace details_spmv{
    template<class Field>
    struct Coo{
        private:
            using Self = Coo<Field>;
        public:            
            typename Field::Element val = 0;
            index_t col = 0;
            index_t row = 0;
            
            Coo() = default;
            Coo(typename Field::Element v, index_t r, index_t c) : val(v), col(c), row(r)
            {}         
            Coo(const Self &) = default;
            Coo(Self &&) = default;
            
            Self& operator=(const Self&) = default;
            Self& operator=(Self &&) = default;
    };
}

// TODO : faster version using fscanf
template <class Field>
void readSmsFormat(const std::string& path, const Field& f,
                   index_t *& row, index_t *& col,
                   typename Field::Element_ptr& val,
                   index_t& rowdim, index_t& coldim, uint64_t & nnz) {
    using namespace details_spmv;
    std::ifstream file(path, std::ios::in);
    std::vector<std::string> tokens;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::copy(std::istream_iterator<std::string>(is),
            std::istream_iterator<std::string>(),
            std::back_inserter<std::vector<std::string> >(tokens));
    rowdim = static_cast<index_t>(std::stoull(tokens[0]));
    coldim = static_cast<index_t>(std::stoull(tokens[1]));
    std::vector<Coo<Field>> data;
    nnz = 0;
    while (std::getline(file, line)) {
        tokens.resize(0);
        std::istringstream iss(line);

        std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string> >(tokens));
    
        if (!(tokens[0] == "0" && tokens[1] == "0" && tokens[2] == "0")) {
            typename Field::Element v; 
            f.init(v, std::stol(tokens[2]));
            index_t r = (index_t)(std::stoull(tokens[0])) - 1;
            index_t c = (index_t)(std::stoull(tokens[1])) - 1;
            data.emplace_back(v, r, c);
        }
    }
    std::sort(data.begin(), data.end(), [](const Coo<Field> & a, const Coo<Field> & b){
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
    });
    row = fflas_new<index_t>(data.size());
    col = fflas_new<index_t>(data.size());
    val = fflas_new(f, data.size(), 1);
    nnz = data.size();
    cout << "nnz : " << nnz << endl;
    for(size_t i = 0, end = data.size() ; i < end ; ++i){
      val[i] = data[i].val;
      col[i] = data[i].col;
      row[i] = data[i].row;
    }
}

int main(int argc, char** argv){
    // using Field = UnparametricField<double>;
    using Field = Givaro::Modular<double>;
    // using Field = class RNSInteger;
    using Element = typename Field::Element;

    
    Field F(11);
    
    std::string path;
       
    COO<Field> Mat;
    
    // mpolyout2.sms
    path = "matrix/EX6.sms";
    if(argc > 1)
        path = argv[1];
       
    readSmsFormat(path, F, Mat.row, Mat.col, Mat.dat, Mat.m, Mat.n, Mat.z);

    std::vector<size_t> tmpv(Mat.m, 0);
    for(size_t i = 0 ; i < Mat.m ; ++i)
        tmpv[Mat.row[i]]++;
    Mat.maxrow = *(std::max_element(tmpv.begin(), tmpv.end()));

    
    VECT<Field> x, y, y2;
    int blockSize = 4;
    cout << "Mat " << Mat.m << " " << Mat.n << endl;// " ; Mat 2 " << Mat2.m << " " << Mat2.n << endl;
    x.dat = fflas_new(F, Mat.n*blockSize, 1, Alignment::CACHE_PAGESIZE);
    y.dat = fflas_new(F, Mat.m*blockSize, 1, Alignment::CACHE_PAGESIZE);


    cout << "fllas_new ok" << endl;
    
    for(size_t i = 0 ; i < Mat.m*blockSize ; ++i)
    {
        x.dat[i] = 1;
    }
    
    for(size_t i = 0 ; i < Mat.n*blockSize ; ++i)
    {
        y.dat[i] = 0;
    }

    fspmm(F, Mat, blockSize, x, blockSize, 1, y, Mat.n);
    cout << "Mat ok" << endl;
    Timer t;
    
    for(size_t i = 0 ; i < 12 ; ++i)
    {
        cout << y.dat[i] << " ";
    }    
    cout << endl;
    
    for(size_t i = 0 ; i < 12 ; ++i)
    {
        
    }    
    cout << endl;
    
    // cout << ((std::equal(y.dat, y.dat+Mat.m, y2.dat)) ? "CORRECT" : "ERROR") << endl;

    sp_delete(Mat);
    // sp_delete(Mat2);
    
    return 0;
}
