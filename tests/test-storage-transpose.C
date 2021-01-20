/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *           
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
#include <iomanip>
#include <iostream>
#include <random>

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/fflas_memory.h"

#include <givaro/modular.h>
#include <recint/rint.h>

#include "fflas-ffpack/fflas/fflas_transpose.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;

using Givaro::Modular;
using Givaro::ModularBalanced;


template <class Field>
bool check_transpose(const Field& F, size_t m, size_t n,
		     const typename Field::Element_ptr A, size_t lda,
		     const typename Field::Element_ptr At, size_t ldat)
{

  // cerr<<"A="<<endl;
  // WriteMatrix (cerr,F,m,n,A,lda);
  // cerr<<"At="<<endl;
  // WriteMatrix (cerr,F,n,m,At,ldat);
  for (size_t i=0;i<m;i++)
    for (size_t j=0;j<n;j++)
      if (A[i*lda+j]!=At[j*ldat+i]) {
	cerr<<"FFLAS::transpose matrix Error at ("<<i<<","<<j<<") -> "<<A[i*lda+j]<<" <> "<<At[j*ldat+i]<<endl;
	cerr<<"A="<<endl;
	WriteMatrix (cerr,F,m,n,A,lda);
	cerr<<"At="<<endl;
	WriteMatrix (cerr,F,n,m,At,ldat);
	return false;
      }
  return true;
}
template <class Field>
bool check_equal(const Field& F, size_t m, size_t n,
		 const typename Field::Element_ptr A, size_t lda,
		 const typename Field::Element_ptr At, size_t ldat)
{


  for (size_t i=0;i<m;i++)
    for (size_t j=0;j<n;j++)
      if (A[i*lda+j]!=At[i*ldat+j]) {
	cerr<<"FFLAS:: equal matrix Error at ("<<i<<","<<j<<") -> "<<A[i*lda+j]<<" <> "<<At[i*ldat+j]<<endl;
	cerr<<"A="<<endl;
	WriteMatrix (cerr,F,m,n,A,lda);
	cerr<<"At="<<endl;
	WriteMatrix (cerr,F,m,n,At,ldat);
	return false;
      }
  return true;
}

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, size_t m, size_t n,  size_t iters, size_t seed, bool bench){
  bool ok = true ;

  int nbit=(int)iters;
  double timeSIMD=0.0, timeNOSIMD=0.0, timeNAIVE=0.0, timeINPLACE=0.0;
  srand(seed);
  while (ok &&  nbit){
    typedef typename Field::Element_ptr Element_ptr;
    // choose Field
    Field* F= chooseField<Field>(q,b,seed);
    if (F==nullptr)
      return true;

    std::ostringstream oss;
    F->write(oss);
    std::cout.fill('.');
    std::cout<<"Checking ";
    std::cout.width(50);
    std::cout<<oss.str();
    std::cout<<" ... ";


    typename Field::RandIter R(*F,seed++);

    /////////////////////////////////////////////////////
    // CHECKING TRANSPOSE (m can be different from n)) //
    /////////////////////////////////////////////////////
    Element_ptr A  = fflas_new(*F,m,n);
    Element_ptr At = fflas_new(*F,m,n);
    Element_ptr Att = fflas_new(*F,m,n);
    RandomMatrix(*F, m, n, A,n, R);

    ftranspose(*F,m,n,A,n,At,m);
    ok&= check_transpose(*F,m,n,A,n,At,m);
    
    ftranspose(*F,n,m,At,m,Att,n); 
    ok&= check_transpose(*F,n,m,At,m,Att,n);
    ok&= check_equal(*F, m, n, A, n, Att, n);
   
    // with submatrix    
    size_t m1 = rand()% m;
    size_t n1 = rand()% n;
    size_t i1 = rand() %(m-m1);
    size_t j1 = rand() %(n-n1);
    /* cout<<"m1="<<m1<<" \nn1="<<n1<<endl; */
    /* cout<<"i1="<<i1<<" \nj1="<<j1<<endl; */
    Element_ptr Abis=A+i1*n+j1;
    Element_ptr Attbis=Att+i1*n+j1;

    ftranspose(*F,m1,n1,Abis,n,At,m1);
    ok&= check_transpose(*F,m1,n1,Abis,n,At,m1);

    ftranspose(*F,n1,m1,At,m1,Attbis,n);
    ok&= check_transpose(*F,n1,m1,At,m1,Attbis,n);
    ok&= check_equal(*F, m1, n1, Abis, n, Attbis,n);
   
    /////////////////////////////////////////////////////
    // CHECKING TRANSPOSEIN (m must equal n))          //
    /////////////////////////////////////////////////////
    
    Element_ptr Amm  = fflas_new(*F,m,m);
    Element_ptr Ammt = fflas_new(*F,m,m);
    RandomMatrix(*F, m, m, Amm, m, R);

    ftranspose(*F,m,m,Amm,m,Ammt,m);
    ftranspose(*F,m,m,Ammt,m,Ammt,m);
    ok&= check_equal(*F, m, m, Amm, m, Ammt,m);

#if 0
    // with submatrix
    m1 = rand()% m;
    i1 = rand() %(m-m1);
    j1 = rand() %(m-m1);
    Abis=Amm+i1*m+j1;
    ftranspose(*F,m1,m1,Abis,m,Ammt,m);
    ftransposein(*F,m1,m1,Ammt,m);
    ok&= check_equal(*F, m1, m1, Abis, m, Ammt,m);


    m1 = rand()% m;
    n1=  rand()% m;
    i1 = rand() %(m-m1);
    j1 = rand() %(m-m1);
    cout<<endl<<"m1="<<m1<<" \nn1="<<n1<<endl;
    cout<<"i1="<<i1<<" \nj1="<<j1<<endl;

    Abis=Amm+i1*m+j1;
    ftranspose(*F,m1,n1,Abis,m,Ammt,m);
    WriteMatrix (cerr,*F,n1,m1,Ammt,m);
    ftransposein(*F,n1,m1,Ammt,m);
    ok&= check_equal(*F, m1, n1, Abis, m, Ammt,m);
#endif
      
    nbit--;
    if ( !ok )
      {std::cout << "FAILED "<<std::endl; return false;}
    else
      std::cout << "PASSED "<<std::endl;

    if (bench){
      Givaro::Timer chrono;
      chrono.clear(); chrono.start();
      ftranspose<Field> (*F,m,n,A,n,At,m);
      chrono.stop();
      timeSIMD+=chrono.usertime();
      chrono.clear(); chrono.start();
      ftranspose<Field, NoSimd<typename Field::Element>> (*F,m,n,A,n,At,m);
      chrono.stop();
      timeNOSIMD+=chrono.usertime();

      /* Naive implem */
      chrono.clear(); chrono.start();
      for(size_t i=0;i<m;i++)
        for(size_t j=0;j<n;j++)
            At[j*m+i]=A[i*n+j];
      chrono.stop();
      timeNAIVE+=chrono.usertime();

      /* inplace */
      chrono.clear(); chrono.start();
      ftranspose<Field> (*F,m,m,Amm,m,Amm,m);
      chrono.stop();
      timeINPLACE+=chrono.usertime();
    }
    delete F;
    fflas_delete(A);
    fflas_delete(At);
    fflas_delete(Att);
    fflas_delete(Amm);
    fflas_delete(Ammt);
  }

  if (bench){
    std::cerr<<std::setprecision(4)<<" running time  naive: "<<timeNAIVE<<std::endl;
    std::cerr<<std::setprecision(4)<<" running time NOSIMD: "<<timeNOSIMD<<std::endl;
    std::cerr<<std::setprecision(4)<<" running time   SIMD: "<<timeSIMD<<std::endl;
    std::cerr<<std::setprecision(4)<<" running time inplace: "<<timeINPLACE<<std::endl;
  }
    
  return ok;   
}



int main(int argc, char** argv)
{
  std::cout<<std::setprecision(17);
  std::cerr<<std::setprecision(17);
  uint64_t seed = getSeed();
  size_t iters = 3 ;
  Givaro::Integer q = -1 ;
  uint64_t b = 0 ;
  size_t m = 50 ;
  size_t n = 50 ;
  bool loop = false;
  bool time= false;
  Argument as[] = {
		   { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		   { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
		   { 'm', "-m M", "Set the row dimension of the matrix.",      TYPE_INT , &m },
		   { 'n', "-n N", "Set the column dimension of the matrix.",      TYPE_INT , &n },
		   { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		   { 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
		   { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
		   { 't', "-t Y/N", "run the test and monitor timing ", TYPE_BOOL, &time },
		   END_OF_ARGUMENTS
  };

  FFLAS::parseArguments(argc,argv,as);
  bool ok = true;
  do{
    ok = ok && run_with_field<Modular<double> >(q,b,m,n,iters, seed, time);
    ok = ok && run_with_field<Modular<float> >(q,b,m,n,iters, seed, time);
    ok = ok && run_with_field<Modular<int64_t> >(q,b,m,n,iters, seed, time);
    ok = ok && run_with_field<Modular<int32_t> >(q,b,m,n,iters, seed, time);
    ok = ok && run_with_field<Modular<int16_t> >(q,b,m,n,iters, seed, time);    
    // ok = ok && run_with_field<Modular<RecInt::rint<7> > >(q,b?b:63_ui64,m,n,iters, seed, time);
    // ok = ok && run_with_field<Modular<RecInt::rint<8> > >(q,b?b:127_ui64,m,n,iters, seed, time);
    // ok = ok && run_with_field<Modular<RecInt::ruint<7>,RecInt::ruint<8> > >(q,b?b:127_ui64,m,n,iters, seed, time);
    /* ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512_ui64),m,n,iters,seed, time); */
    /* ok = ok && run_with_field<Givaro::ZRing<Givaro::Integer> >(0,(b?b:512_ui64),m,n,iters,seed, time); */

    seed++;
  } while (loop && ok);

  return !ok ;

}
