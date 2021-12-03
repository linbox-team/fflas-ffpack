/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Hippolyte Signargout
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


//-------------------------------------------------------------------------
//      Test suite for the Quasi-Separable matrices in SSS format
//-------------------------------------------------------------------------

/* Structure taken from test-quasisep.C */

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

#include <random>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

template<class Field>
bool test_diag_product (const Field & F, size_t n, size_t s, size_t t,
			typename Field::ConstElement_ptr D, size_t ldd,
			typename Field::Element_ptr TS, size_t ldt)
{
  typedef typename Field::Element_ptr Element_ptr ;
  Element_ptr null = fflas_new (F, n, s);
  fzero (F, n, s, null, s);
  Element_ptr C = fflas_new (F, n, t);
    
  productSSSxTS (F, n, t, s, F.one, null, s, null, s, null, s, null, s, null, s, null, s,
		 D, ldd, t, TS, ldt, F.zero, C, t);

  std::cout<<"D:"<<std::endl;
  for (size_t line = 0; line < n; line++)
    {
      for (size_t column = 0; column < s; column++)
	{
	  std::cout<<D[column + line*s]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"B:"<<std::endl;
  for (size_t line = 0; line < n; line++)
    {
      for (size_t column = 0; column < t; column++)
	{
	  std::cout<<TS[column + line*ldt]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"C:"<<std::endl;
  for (size_t line = 0; line < n; line++)
    {
      for (size_t column = 0; column < t; column++)
	{
	  std::cout<<C[column + line*t]<<" ";
	}
      std::cout<<std::endl;
    }
  FFLAS::fflas_delete(null);
  FFLAS::fflas_delete(C);
  return true;
}

/* @brief Tests SSSxTS product on a lower triangular matrix */
template<class Field>
bool test_LT_product (const Field & F, size_t n, size_t s, size_t t,
		      typename Field::ConstElement_ptr P, size_t ldp,
		      typename Field::ConstElement_ptr Q, size_t ldq,
		      typename Field::ConstElement_ptr R, size_t ldr,
		      typename Field::Element_ptr TS, size_t ldt)
{
  typedef typename Field::Element_ptr Element_ptr ;
  Element_ptr null = fflas_new (F, n, s);
  fzero (F, n, s, null, s);
  Element_ptr C = fflas_new (F, n, t);
    
  productSSSxTS (F, n, t, s, F.one, P, ldp, Q, ldq, R, ldr, null, s, null, s, null, s,
		 null, s, TS, ldt, F.zero, C, t);

  std::cout<<"P:"<<std::endl;
  for (size_t line = 0; line < n - s; line++)
    {
      for (size_t column = 0; column < s; column++)
	{
	  std::cout<<P[column + line*ldp]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"Q:"<<std::endl;
  for (size_t line = 0; line < n - s; line++)
    {
      for (size_t column = 0; column < s; column++)
	{
	  std::cout<<Q[column + line*ldq]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"R:"<<std::endl;
  for (size_t line = 0; (line < n - 2*s) && (2*s <= n); line++)
    {
      for (size_t column = 0; column < s; column++)
	{
	  std::cout<<R[column + line*ldq]<<" ";
	}
      std::cout<<std::endl;
    }
    
  std::cout<<"B:"<<std::endl;
  for (size_t line = 0; line < n; line++)
    {
      for (size_t column = 0; column < t; column++)
	{
	  std::cout<<TS[column + line*ldt]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"C:"<<std::endl;
  for (size_t line = 0; line < n; line++)
    {
      for (size_t column = 0; column < t; column++)
	{
	  std::cout<<C[column + line*t]<<" ";
	}
      std::cout<<std::endl;
    }
    
  /* Free memory and return */
  FFLAS::fflas_delete(null);
  FFLAS::fflas_delete(C);
  return true;
}

/* @brief Tests SSSxTS product on an upper triangular matrix */
template<class Field>
bool test_UT_product (const Field & F, size_t n, size_t s, size_t t,
		      typename Field::ConstElement_ptr P, size_t ldp,
		      typename Field::ConstElement_ptr Q, size_t ldq,
		      typename Field::ConstElement_ptr R, size_t ldr,
		      typename Field::Element_ptr TS, size_t ldt)
{
  typedef typename Field::Element_ptr Element_ptr ;
  Element_ptr null = fflas_new (F, n, s);
  fzero (F, n, s, null, s);
  Element_ptr C = fflas_new (F, n, t);

  for (int dec = 0; dec < 4; dec++)
    {
      productSSSxTS (F, n - dec, t, s, F.one, null, s, null, s, null, s, P, ldp, Q, ldq, R, ldr, 
		     null, s, TS, ldt, F.zero, C, t);

      std::cout<<"U:"<<std::endl;
      for (size_t line = 0; line < n - dec- s; line++)
	{
	  for (size_t column = 0; column < s; column++)
	    {
	      std::cout<<P[column + line*ldp]<<" ";
	    }
	  std::cout<<std::endl;
	}
      std::cout<<"V:"<<std::endl;
      for (size_t line = 0; line < n - dec - s; line++)
	{
	  for (size_t column = 0; column < s; column++)
	    {
	      std::cout<<Q[column + line*ldq]<<" ";
	    }
	  std::cout<<std::endl;
	}
      std::cout<<"W:"<<std::endl;
      for (size_t line = 0; (line < n - dec - 2*s) && (2*s <= n - dec); line++)
	{
	  for (size_t column = 0; column < s; column++)
	    {
	      std::cout<<R[column + line*ldq]<<" ";
	    }
	  std::cout<<std::endl;
	}
	
      std::cout<<"B:"<<std::endl;
      for (size_t line = 0; line < n - dec; line++)
	{
	  for (size_t column = 0; column < t; column++)
	    {
	      std::cout<<TS[column + line*ldt]<<" ";
	    }
	  std::cout<<std::endl;
	}
      std::cout<<"C:"<<std::endl;
      for (size_t line = 0; line < n - dec; line++)
	{
	  for (size_t column = 0; column < t; column++)
	    {
	      std::cout<<C[column + line*t]<<" ";
	    }
	  std::cout<<std::endl;
	}
    }
    
  /* Free memory and return */
  FFLAS::fflas_delete(null);
  FFLAS::fflas_delete(C);
  return true;
}

template<class Field>
bool hand_test_reconstruction (const Field & F, size_t n, size_t s, 
			       typename Field::ConstElement_ptr P, size_t ldp,
			       typename Field::ConstElement_ptr Q, size_t ldq,
			       typename Field::ConstElement_ptr R, size_t ldr,
			       typename Field::ConstElement_ptr U, size_t ldu,
			       typename Field::ConstElement_ptr V, size_t ldv,
			       typename Field::ConstElement_ptr W, size_t ldw,
			       typename Field::ConstElement_ptr D, size_t ldd)
{
  typedef typename Field::Element_ptr Element_ptr ;
  Element_ptr C = fflas_new (F, n, n);
   
  for (int dec = 0; dec < 4; dec++)
    {
      sssToDense (F, n - dec, s, P, ldp, Q, ldq, R, ldr, U, ldu, V, ldv, W, ldw,
		  D, ldd, C, n - dec);
      std::cout<<"A:"<<std::endl;
      for (size_t line = 0; line < n - dec; line++)
	{
	  for (size_t column = 0; column < n - dec; column++)
	    {
	      std::cout<<C[column + line*(n - dec)]<<"\t";
	    }
	  std::cout<<std::endl;
	}
    }
  FFLAS::fflas_delete(C);
  return true;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t s, size_t t, size_t iters, uint64_t seed){
  bool ok = true ;
  typedef typename Field::Element_ptr Element_ptr ;
  while (ok && iters)
    {
      Field* F= chooseField<Field>(q,b,seed);
      if (F==nullptr)
	return true;
  
      Element_ptr D = fflas_new (*F, n, s);
      Element_ptr P = fflas_new (*F, n, s);     // Could be n - s but used as U (n - rs) too
      Element_ptr Q = fflas_new (*F, n, s);     // Could be n - rs
      Element_ptr R = fflas_new (*F, n - s, s); // Could be n - s - rs
      Element_ptr B = fflas_new (*F, n, t);
  
      for (size_t line = 0; line < n; line++)
	{
	  for (size_t column = 0; column < s; column++)
	    {
	      D[column + line*s] = line;
	      P[column + line * s] = 1;
	      Q[column + line * s] = 1;
	      if (line < n - s)
		{
		  R[column + line * s] = 2;
		}
	    }
	  for (size_t column = 0; column < t; column++)
	    {
	      B[column + line*t] = column + line;
	    }
	}
      /* Test with only Diagonal blocks */
      //  ok = ok && test_diag_product (*F, n, s, t, D, s, B, t);
      /* Test with only lower-triangular part */
      //  ok = ok && test_LT_product (*F, n, s, t, P, s, Q, s, R, s, B, t);
      /* Test with only lower-triangular part */
      //  ok = ok && test_UT_product (*F, n, s, t, P, s, Q, s, R, s, B, t);
      /* Test reconstruction */
      ok = ok && hand_test_reconstruction (*F, n, s, P, s, Q, s, R, s, P, s, Q, s, R, s, D, s);
      if ( !ok )
	std::cout << "FAILED "<<std::endl;
      else
	std::cout << "PASSED "<<std::endl;
      delete F;

      /* Free memory and return */
      FFLAS::fflas_delete(D);
      FFLAS::fflas_delete(B);
      FFLAS::fflas_delete(P);
      FFLAS::fflas_delete(Q);
      FFLAS::fflas_delete(R);
      iters--;
    }
  
  return ok;
}

int main(int argc, char** argv)
{
  cerr<<setprecision(20);
  Givaro::Integer q=101;
  bool ok=true;
    
    
  /* Very easy */
#if 0
  size_t n = 5;
  size_t s = 1;
  size_t t = 1;
#endif
  /* A bit less easy */
  size_t s = 2;
  size_t t = 2;
  size_t n = 12;

  uint64_t b = 1;

  uint64_t seed = getSeed();
  srand(seed);
    
  ok = ok && run_with_field<Givaro::Modular<int64_t> > (q, b, n, s, t, seed);
  return !ok;
}

/* New main function WIP */
#if 0
int main(int argc, char** argv)
{
  cerr<<setprecision(20);
  Givaro::Integer q=-1;
  size_t b=0;
  size_t n=93;
  size_t s=8;
  size_t t=6;
  int iters=3;
  bool loop=false;
  uint64_t seed = getSeed();

  Argument as[] = {
		   { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
		   { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
		   { 'n', "-n N", "Set the matrix order.", TYPE_INT , &n },
		   { 'r', "-r R", "Set the rank.", TYPE_INT , &r },
		   { 'm', "-m M", "Set the col dim of the TS matrix.", TYPE_INT , &m },
		   { 't', "-t T", "Set the order of quasi-separability.", TYPE_INT , &t },
		   { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
		   { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
		   { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
		   END_OF_ARGUMENTS
  };

  parseArguments(argc,argv,as);

  if (t > n) t = n/2;

  srand(seed);
    
  bool ok=true;
  do{
    ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,b,n,s,t,iters,seed);
    ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,b,n,s,t,iters,seed);
    seed++;
  } while (loop && ok);

  if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;
    
  return !ok;
}
#endif
