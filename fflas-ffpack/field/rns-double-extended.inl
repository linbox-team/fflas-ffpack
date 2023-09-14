/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 *.
 */


#ifndef __FFLASFFPACK_field_rns_double_extended_INL
#define __FFLASFFPACK_field_rns_double_extended_INL
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/fflas/fflas_freduce.h"

namespace FFLAS {
  template<>
  void fscalin(const Givaro::ModularExtended<double> & F, const size_t n, double a,		
               double* X, const size_t incX)
  {
    for(size_t i=0;i<n;i+=incX){
      F.mulin(X[i],a);
    }
  }
}

namespace FFPACK {

  /**************************************
   * V1 Variant of RNS DOUBLE EXTENDED
   **************************************/


  // reduce entries of Arns to be less than the rns basis elements
  inline void rns_double_extended_V1::reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR) const{

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    using simd = Simd<double>;
    using vect_t = typename simd::vect_t;

    if(_size % simd::vect_size == 0){
      //#pragma omp parallel for schedule(static, 256)			  
      for(size_t i = 0 ; i < n ; i++){
        vect_t tmp1, tmp2, v, min, max, basis, inv, neg;
        min = simd::set1(0.);
        for(size_t j = 0 ; j < _size ; j+=simd::vect_size){
          basis = simd::load(_basis.data()+j);
          inv   = simd::load(_invbasis.data()+j);
          max   = simd::load(_basisMax.data()+j);
          neg   = simd::load(_negbasis.data()+j);
          v     = simd::load(Arns+i*_size+j);
          simd::mod(v, basis, inv, neg, min, max, tmp1,tmp2);
          simd::store(Arns+i*_size+j, v);
        }
      }
    } else{
      //#pragma omp parallel for schedule(static, 256)			  
      for(size_t i = 0 ; i < n ; i++){
        vect_t tmp1, tmp2, tmp3, v, min, max, basis, inv, neg;
        size_t j = 0;
        for( ; j < ROUND_DOWN(_size, simd::vect_size) ; j+=simd::vect_size){
          basis = simd::load(_basis.data()+j);
          inv   = simd::load(_invbasis.data()+j);
          max   = simd::load(_basisMax.data()+j);
          neg   = simd::load(_negbasis.data()+j);
          v     = simd::load(Arns+i*_size+j);
          simd::mod(v, basis, inv, neg, min, max, tmp1,tmp2);
          simd::store(Arns+i*_size+j, v);
        }
        for( ; j < _size ; ++j){
          _field_rns[j].reduce(Arns[i*_size+j]);
        }
      }
    }
#else
    for(size_t i=0;i<_size;i++)
      FFLAS::freduce (_field_rns[i],n,Arns+i*rda,1);

#endif

  }



  void rns_double_extended_V1::init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR) const
  {
    if (_size>= (1<<16)){
      std::cerr<<"RNS EXTENDED DOUBLE: init Error -> the nbr of moduli in RNS basis is > 2^16, not implemented. aborting\n";std::terminate();
    }
#ifdef BENCH_RNS
    if (m!=1 && n!=1){
      std::cerr<<RALIGN<<"RNS double ext (To) --> rns size ("<<_size<<") kronecker size ("<<k<<") data dim ("<<m*n<<")"<<std::endl;
      std::cerr<<"RNS double ext -> Numbit(M)="<<_M.bitsize()<<std::endl;
    }
#endif
    //init(m*n,Arns,A,lda);
    if (k>_ldm){
      FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_struct: init (too large entry)");
      std::cerr<<"k="<<k<<" _ldm="<<_ldm<<std::endl;
    }
    size_t mn=m*n;
    size_t mnk=mn*k;
    double *A_beta = FFLAS::fflas_new<double >(mnk*6);
    double *A_beta0=A_beta;
    double *A_beta1=A_beta+1*mnk;
    double *A_beta2=A_beta+2*mnk;
    double *A_beta3=A_beta+3*mnk;
    double *A_beta4=A_beta+4*mnk;
    double *A_beta5=A_beta+5*mnk;
    size_t mnsize=mn*_size;
    double *A_rns_tmp = FFLAS::fflas_new<double >(mnsize*6);
    double *A_rns0=A_rns_tmp;
    double *A_rns1=A_rns_tmp+1*mnsize;
    double *A_rns2=A_rns_tmp+2*mnsize;
    double *A_rns3=A_rns_tmp+3*mnsize;
    double *A_rns4=A_rns_tmp+4*mnsize;
    double *A_rns5=A_rns_tmp+5*mnsize;

		
    const integer* Aiter=A;
    // split A into A_beta according to a Kronecker transform in base 2^48
    Givaro::Timer tkr; tkr.start();
    PARFOR1D(i,m,SPLITTER(NUM_THREADS),
             for(size_t j=0;j<n;j++){
               size_t idx=j+i*n;
               const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
               const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
               size_t l=0,h=0;
               // size in base 2^16
               size_t k16=((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2);
               size_t maxs=std::min(k,((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3);// to ensure 32 bits portability
               double a0,a1,a2;
               if (m0[0]->_mp_size >= 0)
                 for (;l<maxs;l++){
                   a0= m0_ptr[h++];//std::cout<<"a0="<<(uint64_t)a0<<std::endl;
                   a1= m0_ptr[h++];//std::cout<<"a1="<<(uint64_t)a1<<std::endl;
                   a2= m0_ptr[h++];//std::cout<<"a2="<<(uint64_t)a2<<std::endl;
                   A_beta0[l+idx*k]= a0;
                   A_beta1[l+idx*k]= a1;
                   A_beta2[l+idx*k]= a2;
                   A_beta3[l+idx*k]= a0+a1;
                   A_beta4[l+idx*k]= a1+a2;
                   A_beta5[l+idx*k]= a0+a1+a2;							  							  
                 }
               else
                 for (;l<maxs;l++){
                   a0= -double(m0_ptr[h++]);
                   a1= -double(m0_ptr[h++]);
                   a2= -double(m0_ptr[h++]);
                   A_beta0[l+idx*k]= a0;
                   A_beta1[l+idx*k]= a1;
                   A_beta2[l+idx*k]= a2;
                   A_beta3[l+idx*k]= a0+a1;
                   A_beta4[l+idx*k]= a1+a2;
                   A_beta5[l+idx*k]= a0+a1+a2;							  							  
                 }
               for (;l<k;l++){
                 a0= (h<k16)?m0_ptr[h++]:0.;
                 a1= (h<k16)?m0_ptr[h++]:0.;
                 a2= (h<k16)?m0_ptr[h++]:0.;
                 A_beta0[l+idx*k]= a0;
                 A_beta1[l+idx*k]= a1;
                 A_beta2[l+idx*k]= a2;
                 A_beta3[l+idx*k]= a0+a1;
                 A_beta4[l+idx*k]= a1+a2;
                 A_beta5[l+idx*k]= a0+a1+a2;							  							  							  
               }					
             }
             );

#ifdef CHECK_RNS
    for (size_t i=0;i<m;i++)
      for (size_t j=0;j<n;j++){
        //std::cout<<"A="<<A[i*lda+j]<<std::endl;
				
				
        integer tmp=0;
        int idx=j+i*n;
        for(int64_t l=(int64_t)k-1;l>=0;l--){
          int64_t limb,c0,c1,c2,c3,c4;
          c0= A_beta0[l+idx*k]; // 1
          c1= A_beta3[l+idx*k] - A_beta1[l+idx*k] - A_beta0[l+idx*k] ; // X
          c2= A_beta5[l+idx*k] - A_beta3[l+idx*k] - A_beta4[l+idx*k] + 2*A_beta1[l+idx*k]; // X^2
          c3= A_beta4[l+idx*k] - A_beta1[l+idx*k] - A_beta2[l+idx*k]; // X^3
          c4= A_beta2[l+idx*k]; // X^4
          if (c1!=0. || c3!=0.) std::cout<<"RNS EXTENDED: error in splitting entries (linear form)\n";
          limb= c0+(c2<<16)+(c4<<32);
          tmp=(tmp<<48)+limb;
        }
        if (tmp!=A[i*lda+j]){
          std::cout<<"RNS EXTENDED: error in splitting entries (16-adic) --> "<<tmp<<" != "<<A[i*lda+j]<<std::endl;;
          std::cout<<"MAX="<<std::min(k,((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3)<<std::endl;
          std::cout<<"MAX1="<<k<<std::endl;
          std::cout<<"MAX2="<<((Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2)/3<<std::endl;
          std::cout<<"bisize: "<<Aiter[j+i*lda].bitsize()<<std::endl;

          tmp=A[i*lda+j];
          for(size_t l=0;l<k;l++){
            std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFF))             << " != "<<int64_t(A_beta0[l+idx*k])<<std::endl;
            std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFFFFFF)>>16)     << " != "<<int64_t(A_beta1[l+idx*k])<<std::endl;
            std::cout<<"l="<<l<<"  -->  "<<((tmp[0] & 0xFFFFFFFFFFFF)>>32) << " != "<<int64_t(A_beta2[l+idx*k])<<std::endl;
            tmp>>=48;
          }
        }
      }													
#endif
    tkr.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<RALIGN<<"RNS double ext (To) - Kronecker : "<<tkr.usertime()<<std::endl;
#endif


    Givaro::ZRing<double> ZD;

    // Using Helper for potential parallelism -> need to be activated by hand
    FFLAS::MMHelper<Givaro::ZRing<double>, FFLAS::MMHelperAlgo::Winograd>  MMH (ZD, -1, FFLAS::ParSeqHelper::Sequential());

    if (RNS_MAJOR==false) {
      // A_rns = _crt_in x A_beta^T
      Givaro::Timer tfgemm; tfgemm.start();
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[0].data(),_ldm,A_beta0,k,0.,A_rns0, mn, MMH);
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[1].data(),_ldm,A_beta1,k,0.,A_rns1, mn, MMH);
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[2].data(),_ldm,A_beta2,k,0.,A_rns2, mn, MMH);
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[3].data(),_ldm,A_beta3,k,0.,A_rns3, mn, MMH);
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[4].data(),_ldm,A_beta4,k,0.,A_rns4, mn, MMH);
      FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in[5].data(),_ldm,A_beta5,k,0.,A_rns5, mn, MMH);

      tfgemm.stop();
#ifdef BENCH_RNS
      if (m!=1 && n!=1)
        std::cerr<<RALIGN<<"RNS double ext (To)  - fgemm : "<<tfgemm.usertime()<<std::endl;
#endif
    }
    else {
      // Arns =  A_beta x _crt_in^T
      std::cerr<<"NOT YET IMPLEMENTED .... aborting\n"; std::terminate();
    }
    Givaro::Timer tred; tred.start();

#ifdef RNS_DEBUG			
    std::cout<<"A:=Matrix("<<m<<","<<n<<",[";
    for(size_t i=0;i<m;i++){
      std::cout<<"[";
      for(size_t j=0;j<n;j++)
        std::cout<<(int64_t)A[j+i*lda]<<(j!=n-1?",":(i!=m-1?"],":"]"));
    }
    std::cout<<"]);\n";			

    for(size_t l=0;l<6;l++){
      std::cout<<"Abeta"<<l<<":=Matrix("<<_size<<","<<mn<<",[";
      for(size_t i=0;i<k;i++){
        std::cout<<"[";
        for(size_t j=0;j<mn;j++)
          std::cout<<(int64_t)A_beta[l*mnk+j+i*mn]<<(j!=mn-1?",":(i!=k-1?"],":"]"));
      }
      std::cout<<"]);\n";
    }

    for(size_t l=0;l<6;l++){
      std::cout<<"Arns"<<l<<":=Matrix("<<_size<<","<<mn<<",[";
      for(size_t i=0;i<_size;i++){
        std::cout<<"[";
        for(size_t j=0;j<mn;j++)
          std::cout<<(int64_t)A_rns_tmp[l*mnsize+j+i*mn]<<(j!=mn-1?",":(i!=_size-1?"],":"]"));
      }
      std::cout<<"]);\n";
    }
#endif
			
			
    double c0,c1,c2,c3,c4;
    //double C,C1,C3;
    //int64_t c0,c1,c2,c3,c4;

    const double two16=1UL<<16,
      two32=1UL<<32,
      two48=1UL<<48,
      two64=double(1UL<<32)*double(1UL<<32); 
			
    PARFOR1D(i,_size,SPLITTER(NUM_THREADS),
             double two16_mod_mi,two48_mod_mi;
             _field_rns[i].init(two16_mod_mi,(1<<16));
             _field_rns[i].mul(two48_mod_mi,two16_mod_mi,two16_mod_mi);
             _field_rns[i].mulin(two48_mod_mi,two16_mod_mi);
					 
             double two_16_invmi = two16 * _invbasis[i];
             double two_32_invmi = two32 * _invbasis[i];
             double two_48_invmi = two48 * _invbasis[i]; 
             double two_64_invmi = two64 * _invbasis[i];
             double q1,q2,q3,q4;

             // std::cout<<"M1:="<<two_16_invmi<<";\n";
             // std::cout<<"M2:="<<two_32_invmi<<";\n";
             // std::cout<<"M3:="<<two_48_invmi<<";\n";
             // std::cout<<"M4:="<<two_64_invmi<<";\n";
					 
					 
             for(size_t j=0;j<mn;j++)
               {
                 size_t h=j+i*mn;
                 c0= A_rns0[h]; // 1
                 c1= A_rns3[h] - A_rns1[h] - A_rns0[h] ; // X
                 c2= A_rns5[h] - A_rns3[h] - A_rns4[h] + 2*A_rns1[h]; // X^2
                 c3= A_rns4[h] - A_rns1[h] - A_rns2[h]; // X^3
                 c4= A_rns2[h]; // X^4
#ifdef RNS_DEBUG
                 std::cout<<(int64_t)c0
                          <<"+2^16*"<<(int64_t)c1
                          <<"+2^32*"<<(int64_t)c2
                          <<"+2^48*"<<(int64_t)c3
                          <<"+2^64*"<<(int64_t)c4<<";\n";
#endif				 
                 // compute c=c0+c1.2^16+c2.2^32+c3.2^48+c4.2^64 mod mi							 
                 // _field_rns[i].axpy(c,c4,two16_mod_mi,c3);
                 // _field_rns[i].axpy(c,c ,two16_mod_mi,c2);
                 // _field_rns[i].axpy(c,c ,two16_mod_mi,c1);
                 // _field_rns[i].axpy(c,c ,two16_mod_mi,c0);
                 // Arns[j+i*rda]= c+(c>=0?0:_basis[i]);

                 q1= std::floor(c1*two_16_invmi);//std::cout<<"Q1:="<<(int64_t)q1<<";\n";
                 q2= std::floor(c2*two_32_invmi);//std::cout<<"Q2:="<<(int64_t)q2<<";\n";
                 q3= std::floor(c3*two_48_invmi);//std::cout<<"Q3:="<<(int64_t)q3<<";\n";
                 q4= std::floor(c4*two_64_invmi);//std::cout<<"Q4:="<<(int64_t)q4<<";\n";

                 c1*=two16;
                 c2*=two32;
                 c3*=two48;
                 c4*=two64;
							 
                 c1=fma(q1,_negbasis[i],c1); //std::cout<<"D1:="<<(int64_t)c1<<";\n";
                 c2=fma(q2,_negbasis[i],c2); //std::cout<<"D2:="<<(int64_t)c2<<";\n";
                 c3=fma(q3,_negbasis[i],c3); //std::cout<<"D3:="<<(int64_t)c3<<";\n";
                 c4=fma(q4,_negbasis[i],c4); //std::cout<<"D4:="<<(int64_t)c4<<";\n";

                 c0+=c1;
                 c2+=c3;
                 c0+=c4+c2;
                 while(c0<0.)        c0+=_basis[i];
                 while(c0>_basis[i]) c0-=_basis[i];
                 Arns[j+i*rda]= c0;
							 
                 // c1+=c2<<16;
                 // c3+=c4<<16;
                 // _field_rns[i].reduce(C1,c1);
                 // _field_rns[i].reduce(C3,c3);
                 // _field_rns[i].axpy(C1,C1,two16_mod_mi,c0);
                 // _field_rns[i].mul(C3,C3,two48_mod_mi);
                 // _field_rns[i].add(C,C1,C3);
                 // Arns[j+i*rda]= C+(C>=0?0:_basis[i]);
               }
             );
    //reduce(mn,Arns,rda,RNS_MAJOR);
			
    tred.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<RALIGN<<"RNS double ext (To)  - Reduce : "<<tred.usertime()<<std::endl;
#endif
	
    FFLAS::fflas_delete(A_beta);
    FFLAS::fflas_delete(A_rns_tmp);

#ifdef CHECK_RNS
    bool ok=true;
    for (size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
        for(size_t k=0;k<_size;k++){
          ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0)) == (int64_t) Arns[i*n+j+k*rda]);
          if (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
              != (int64_t) Arns[i*n+j+k*rda])
            {
              std::cout<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
                       <<" != "
                       <<(int64_t) Arns[i*n+j+k*rda]
                       <<" --> "<<A[i*lda+j]<<"  mod ("<<(uint64_t)_basis[k]<<")"<<std::endl;
            }
          else{
            // std::cout<<"OK -> "<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
            // 		 <<" == "
            // 		 <<(int64_t) Arns[i*n+j+k*rda]
            // 		 <<" --> "<<A[i*lda+j]<<"  mod ("<<(uint64_t)_basis[k]<<")"<<std::endl;
						
          }
        }
    std::cout<<"RNS EXTENDED Double (To)  ... "<<(ok?"OK":"ERROR")<<std::endl;
    if (!ok) std::terminate();
#endif



			
  }




  struct uint48_t {
    uint64_t data :48;
    uint48_t () =default;
    uint48_t (uint64_t x) :data(x){}
  } __attribute__((packed));
		

		
  void rns_double_extended_V1::convert(size_t m, size_t n, integer gamma, integer* A, size_t lda, const double* Arns, size_t rda, bool RNS_MAJOR) const {
    // assume the number of moduli is less than 2^16
    if (_size>= (1<<16)){
      std::cerr<<"RNS EXTENDED DOUBLE: convert Error -> the nbr of moduli in RNS basis is > 2^16, not implemented. aborting\n";std::terminate();
    }
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<RALIGN<<"RNS double ext (From) --> rns size ("<<_size<<") kronecker size ("<<_ldm<<") data dim ("<<m*n<<")"<<std::endl;
#endif
			
#ifdef CHECK_RNS
    integer* Acopy=new integer[m*n];
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
        Acopy[i*n+j]=A[i*lda+j];			
#endif

    integer hM= (_M-1)>>1;
    size_t  mn= m*n;
    size_t mnldm=mn*_ldm;
    double *A_beta= FFLAS::fflas_new<double>(6*mnldm);			
    double *A_beta0=A_beta;
    double *A_beta1=A_beta+1*mnldm;
    double *A_beta2=A_beta+2*mnldm;
    double *A_beta3=A_beta+3*mnldm;
    double *A_beta4=A_beta+4*mnldm;
    double *A_beta5=A_beta+5*mnldm;
    size_t mnsize=mn*_size;
    double *A_rns_tmp = FFLAS::fflas_new<double >(mnsize*6);
    double *A_rns0=A_rns_tmp;
    double *A_rns1=A_rns_tmp+1*mnsize;
    double *A_rns2=A_rns_tmp+2*mnsize;
    double *A_rns3=A_rns_tmp+3*mnsize;
    double *A_rns4=A_rns_tmp+4*mnsize;
    double *A_rns5=A_rns_tmp+5*mnsize;
			
    Givaro::Timer tsplit;
    tsplit.start();			
    uint64_t tmp,idx;
    double aa0,aa1,aa2;
    for (size_t i=0;i<_size;i++)
      for(size_t j=0;j<mn;j++){
        idx=i*mn+j;
        tmp=Arns[idx];
        aa0=tmp&0xFFFF;
        aa1=(tmp>>16)&0xFFFF;
        aa2=tmp>>32;
        A_rns0[idx]=aa0;
        A_rns1[idx]=aa1;
        A_rns2[idx]=aa2;
        A_rns3[idx]=aa0+aa1;
        A_rns4[idx]=aa1+aa2;
        A_rns5[idx]=aa0+aa1+aa2;
      }
    tsplit.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)				
      std::cerr<<RALIGN<<"RNS double ext (From) -  split : "<<tsplit.usertime()<<std::endl;
#endif
    Givaro::Timer tfgemmc;tfgemmc.start();
    if (RNS_MAJOR==false){
      // compute A_beta = Ap^T x M_beta
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns0,(int) rda, _crt_out[0].data(),(int) _ldm, 0., A_beta0,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns1,(int) rda, _crt_out[1].data(),(int) _ldm, 0., A_beta1,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns2,(int) rda, _crt_out[2].data(),(int) _ldm, 0., A_beta2,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns3,(int) rda, _crt_out[3].data(),(int) _ldm, 0., A_beta3,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns4,(int) rda, _crt_out[4].data(),(int) _ldm, 0., A_beta4,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
      FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , A_rns5,(int) rda, _crt_out[5].data(),(int) _ldm, 0., A_beta5,(int)_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
    }
    else {
      // compute A_beta = Ap x M_Beta 
      std::cerr<<"NOT YET IMPLEMENTED .... aborting\n"; std::terminate();
    }
    tfgemmc.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<RALIGN<<"RNS double ext (From) -  fgemm : "<<tfgemmc.usertime()<<std::endl;
#endif
    // compute A using inverse Kronecker transform of A_beta expressed in base 2^48

    FFLAS::fflas_delete( A_rns_tmp);
#ifdef RNS_DEBUG			
    std::cout<<"basis:= [";
    for(size_t i=0;i<_size;i++)
      std::cout<<(int64_t)_basis[i]<<(i!=_size-1?",":"];\n");

    std::cout<<"Arns:=Matrix("<<_size<<","<<mn<<",[";
    for(size_t i=0;i<_size;i++){
      std::cout<<"[";
      for(size_t j=0;j<mn;j++)
        std::cout<<(int64_t)Arns[j+i*mn]<<(j!=mn-1?",":(i!=_size-1?"],":"]"));
    }
    std::cout<<"]);\n";
			
    for(size_t l=0;l<6;l++){
      std::cout<<"Arns"<<l<<":=Matrix("<<_size<<","<<mn<<",[";
      for(size_t i=0;i<_size;i++){
        std::cout<<"[";
        for(size_t j=0;j<mn;j++)
          std::cout<<(int64_t)A_rns_tmp[l*mnsize+j+i*mn]<<(j!=mn-1?",":(i!=_size-1?"],":"]"));
      }
      std::cout<<"]);\n";
    }
			
    for(size_t l=0;l<6;l++){
      std::cout<<"Abeta"<<l<<":=Matrix("<<mn<<","<<_size<<",[";
      for(size_t i=0;i<mn;i++){
        std::cout<<"[";
        for(size_t j=0;j<_ldm;j++)
          std::cout<<(int64_t)A_beta[l*mnldm+j+i*_ldm]<<(j!=_ldm-1?",":(i!=mn-1?"],":"]"));
      }
      std::cout<<"]);\n";
    }

#endif

    Givaro::Timer tkroc;
    tkroc.start();
			
    integer* Aiter= A;
    size_t k=_ldm;
    size_t k64= (k*48+64)/64+1;
#ifdef RNS_DEBUG
    std::cout<<"kron base 48: "<<k<<std::endl;
    std::cout<<"kron base 64: "<<k64<<std::endl;
#endif
    std::vector<uint16_t> A0_tmp(k64<<2,0),A1_tmp(k64<<2,0),A2_tmp(k64<<2,0),A3_tmp(k64<<2,0),A4_tmp(k64<<2,0);
    uint48_t *A0,*A1,*A2,*A3,*A4;
    A0= reinterpret_cast<uint48_t*>(A0_tmp.data());
    A1= reinterpret_cast<uint48_t*>(A1_tmp.data()+1);
    A2= reinterpret_cast<uint48_t*>(A2_tmp.data()+2);
    A3= reinterpret_cast<uint48_t*>(A3_tmp.data()+3);
    A4= reinterpret_cast<uint48_t*>(A4_tmp.data()+4);
    integer a0,a1,a2,a3,a4,res;
    mpz_t *m0,*m1,*m2,*m3,*m4;
    m0= reinterpret_cast<mpz_t*>(&a0);
    m1= reinterpret_cast<mpz_t*>(&a1);
    m2= reinterpret_cast<mpz_t*>(&a2);
    m3= reinterpret_cast<mpz_t*>(&a3);
    m4= reinterpret_cast<mpz_t*>(&a4);
    mp_limb_t *m0_d,*m1_d,*m2_d,*m3_d,*m4_d;
    m0_d = m0[0]->_mp_d;
    m1_d = m1[0]->_mp_d;
    m2_d = m2[0]->_mp_d;
    m3_d = m3[0]->_mp_d;
    m4_d = m4[0]->_mp_d;
    m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = m4[0]->_mp_alloc = (int) (k64*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
    m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size  = m3[0]->_mp_size  = m4[0]->_mp_size  = (int) (k64*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
    //		auto sp=SPLITTER();
    //		PARFOR1D(i,m,sp,
#ifdef RNS_DEBUG
    std::cout<<"M:="<<_M<<std::endl;
#endif
    uint64_t c0,c1,c2,c3,c4;
    for(size_t i=0;i<m;i++)
      for (size_t j=0;j<n;j++){
        size_t idx=i*n+j;
#ifdef RNS_DEBUG
        std::cout<<"c0:=0;";
        std::cout<<"c1:=0;";
        std::cout<<"c2:=0;";
        std::cout<<"c3:=0;";
        std::cout<<"c4:=0;";
#endif	
        for (size_t l=0;l<k;l++){
          size_t idxl=l+idx*k;
          c0= A_beta0[idxl]; // 1
          c1= A_beta3[idxl] - A_beta1[idxl] - A_beta0[idxl] ; // X
          c2= A_beta5[idxl] - A_beta3[idxl] - A_beta4[idxl] + 2*A_beta1[idxl]; // X^2
          c3= A_beta4[idxl] - A_beta1[idxl] - A_beta2[idxl]; // X^3
          c4= A_beta2[idxl]; // X^4
#ifdef RNS_DEBUG
          std::cout<<"CCC2"<<l<<":="<<c2<<"\n;";
						
          std::cout<<"c0:=c0+2^"<<l*48<<"*"<<c0<<";";
          std::cout<<"c1:=c1+2^"<<l*48+16<<"*"<<c1<<";";
          std::cout<<"c2:=c2+2^"<<l*48+32<<"*"<<c2<<";";
          std::cout<<"c3:=c3+2^"<<l*48+48<<"*"<<c3<<";";
          std::cout<<"c4:=c4+2^"<<l*48+64<<"*"<<c4<<";";
#endif
          A0[l]= c0;
          A1[l]= c1;
          A2[l]= c2;
          A3[l]= c3;
          A4[l]= c4;
        }
#ifdef RNS_DEBUG
        std::cout<<"0 mod "<<_M<<";\n";
#endif
        // see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
        m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(A0_tmp.data());
        m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(A1_tmp.data());
        m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(A2_tmp.data());
        m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(A3_tmp.data());
        m4[0]->_mp_d= reinterpret_cast<mp_limb_t*>(A4_tmp.data());
        res = a0;res+= a1;res+= a2;res+= a3;res+=a4;
        res%=_M;

#ifdef RNS_DEBUG
        std::cout<<"a0:="<<a0<<";\n";
        std::cout<<"a1:="<<a1<<";\n";
        std::cout<<"a2:="<<a2<<";\n";
        std::cout<<"a3:="<<a3<<";\n";
        std::cout<<"a4:="<<a4<<";\n";
        std::cout<<"res:="<<res<<";\n";
#endif	
        // get the correct result according to the expected sign of A
        if (res>hM)
          res-=_M;
        if (gamma==0)
          Aiter[j+i*lda]=res;
        else
          if (gamma==integer(1))
            Aiter[j+i*lda]+=res;
          else
            if (gamma==integer(-1))
              Aiter[j+i*lda]=res-Aiter[j+i*lda];
            else{
              Aiter[j+i*lda]*=gamma;
              Aiter[j+i*lda]+=res;
            }
					
      }
    tkroc.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<RALIGN<<"RNS double ext (From) - Kronecker : "<<tkroc.usertime()<<std::endl;
#endif
			
    m0[0]->_mp_d = m0_d;
    m1[0]->_mp_d = m1_d;
    m2[0]->_mp_d = m2_d;
    m3[0]->_mp_d = m3_d;
    m4[0]->_mp_d = m4_d;
    m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc= m3[0]->_mp_alloc = m4[0]->_mp_alloc = 1;
    m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size = m3[0]->_mp_size  = m4[0]->_mp_size  = 0;
    FFLAS::fflas_delete( A_beta);
				 
#ifdef CHECK_RNS
    bool ok=true;
    for (size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
        for(size_t k=0;k<_size;k++){
          int64_t _p =(int64_t) _basis[k];
          integer curr=A[i*lda+j] - gamma*Acopy[i*n+j];
          ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);
          if (( curr% _p +(curr%_p<0?_p:0)) != (int64_t) Arns[i*n+j+k*rda])
            std::cout<<"A("<<i<<","<<j<<") -> "<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"!="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
        }
    std::cout<<"RNS EXTENDED Double (From)  ... "<<(ok?"OK":"ERROR")<<std::endl;
    if (!ok) std::terminate();
#endif
  }



  /**************************************
   * V2 Variant of RNS DOUBLE EXTENDED
   **************************************/

  // reduce entries of Arns to be less than the rns basis elements
  inline void rns_double_extended_V2::reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR) const{

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    using simd = Simd<double>;
    using vect_t = typename simd::vect_t;

    if(_size % simd::vect_size == 0){
      //#pragma omp parallel for schedule(static, 256)			  
      for(size_t i = 0 ; i < n ; i++){
        vect_t tmp1, tmp2, v, min, max, basis, inv, neg;
        min = simd::set1(0.);
        for(size_t j = 0 ; j < _size ; j+=simd::vect_size){
          basis = simd::load(_basis.data()+j);
          inv   = simd::load(_invbasis.data()+j);
          max   = simd::load(_basisMax.data()+j);
          neg   = simd::load(_negbasis.data()+j);
          v     = simd::load(Arns+i*_size+j);
          simd::mod(v, basis, inv, neg, min, max, tmp1,tmp2);
          simd::store(Arns+i*_size+j, v);
        }
      }
    } else{
      //#pragma omp parallel for schedule(static, 256)			  
      for(size_t i = 0 ; i < n ; i++){
        vect_t tmp1, tmp2, tmp3, v, min, max, basis, inv, neg;
        size_t j = 0;
        for( ; j < ROUND_DOWN(_size, simd::vect_size) ; j+=simd::vect_size){
          basis = simd::load(_basis.data()+j);
          inv   = simd::load(_invbasis.data()+j);
          max   = simd::load(_basisMax.data()+j);
          neg   = simd::load(_negbasis.data()+j);
          v     = simd::load(Arns+i*_size+j);
          simd::mod(v, basis, inv, neg, min, max, tmp1,tmp2);
          simd::store(Arns+i*_size+j, v);
        }
        for( ; j < _size ; ++j){
          _field_rns[j].reduce(Arns[i*_size+j]);
        }
      }
    }
#else
    for(size_t i=0;i<_size;i++)
      FFLAS::freduce (_field_rns[i],n,Arns+i*rda,1);

#endif

  }


  inline void rns_double_extended_V2::init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR) const
  {
    if (k>_ldm){
      FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_struct: init (too large entry)");
      std::cerr<<"k="<<k<<" _ldm="<<_ldm<<std::endl;
    }
    size_t mn=m*n;
    double *A_beta = FFLAS::fflas_new<double >(mn*k);
    double *Arns_tmp = FFLAS::fflas_new<double >(mn*_size*2);
    const integer* Aiter=A;
    // split A into A_beta according to a Kronecker transform in base 2^16
    //		auto sp=SPLITTER(MAX_THREADS,FFLAS::CuttingStrategy::Column,FFLAS::StrategyParameter::Threads);
    
#ifdef BENCH_RNS
    if (m!=1 && n!=1){
      std::cerr<<"RNS EXTENDED double (To) --> rns size ("<<_size<<") kronecker size ("<<k<<") data dim ("<<m*n<<")"<<std::endl;
      std::cerr<<"RNS EXTENDED double  -> Numbit(M)="<<_M.bitsize()<<std::endl;
    }
#endif    
    Givaro::Timer tkr; tkr.start();
    PARFOR1D(i,m,SPLITTER(NUM_THREADS),
             for(size_t j=0;j<n;j++){
               size_t idx=j+i*n;
               const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
               const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
               size_t l=0;
               //size_t maxs=std::min(k,(Aiter[j+i*lda].size())<<2);
               size_t maxs=std::min(k,(Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2);// to ensure 32 bits portability

               if (m0[0]->_mp_size >= 0)
                 for (;l<maxs;l++)
                   A_beta[l+idx*k]=  m0_ptr[l];
               else
                 for (;l<maxs;l++)
                   A_beta[l+idx*k]= - double(m0_ptr[l]);
               for (;l<k;l++)
                 A_beta[l+idx*k]=  0.;

               // 	   );
             }
             );

    tkr.stop();
#ifdef BENCH_RNS
    if (m!=1 && n!=1)
      std::cerr<<"RNS EXTENDED double (To) - Kronecker : "<<tkr.usertime()<<std::endl;
#endif
    if (RNS_MAJOR==false) {
      // Arns_tmp = _crt_in x A_beta^T
      Givaro::Timer tfgemm; tfgemm.start();
      FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,2*_size,mn,k,1.0,_crt_in.data(),_ldm,A_beta,k,0.,Arns_tmp,mn,
                    FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
      tfgemm.stop();
#ifdef BENCH_RNS
      if(m>1 && n>1) 	std::cerr<<"RNS EXTENDED double (To) - fgemm : "<<tfgemm.usertime()<<std::endl;
#endif
    }
    else {
      // Arns =  A_beta x _crt_in^T
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)mn,(int)_size*2,(int)k,1.0,A_beta,(int)k,_crt_in.data(),(int)_ldm,0.,Arns_tmp,(int)_size*2);
    }
    Givaro::Timer tred; tred.start();
    for(size_t i=0;i<_size;i++){
      FFLAS::fscalin (_field_rns[i],mn, double(1<<_shift), Arns_tmp+(i+_size)*mn,1);
      FFLAS::fadd(_field_rns[i],mn,Arns_tmp+(i)*mn,1,Arns_tmp+(_size+i)*mn,1, Arns+i*rda,1);
      FFLAS::freduce (_field_rns[i],mn, Arns+i*rda,1);
    }
    tred.stop();
#ifdef BENCH_RNS
    if(m>1 && n>1) 			std::cerr<<"RNS EXTENDED double (To) - Reduce : "<<tred.usertime()<<std::endl;
#endif
    FFLAS::fflas_delete( A_beta);
    FFLAS::fflas_delete( Arns_tmp);

#ifdef CHECK_RNS
    bool ok=true;
    for (size_t i=0;i<m;i++)
      for(size_t j=0;j<n;j++)
        for(size_t k=0;k<_size;k++){
          ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0)) == (int64_t) Arns[i*n+j+k*rda]);
          if (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
              != (int64_t) Arns[i*n+j+k*rda])
            {
              std::cout<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
                       <<" != "
                       <<(int64_t) Arns[i*n+j+k*rda]
                       << "    -> "<<A[i*lda+j]
                       <<std::endl;
            }
        }
    std::cout<<"RNS EXTENDED double (To) ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif
  }




  inline void rns_double_extended_V2::convert(size_t m, size_t n, integer gamma, integer* A, size_t lda,
                                              const double* Arns, size_t rda, bool RNS_MAJOR) const
	{
#ifdef BENCH_RNS
	  if (m!=1 && n!=1)
	    std::cerr<<"RNS double ext (From) --> rns size ("<<_size<<") kronecker size ("<<_ldm<<") data dim ("<<m*n<<")"<<std::endl;
#endif

#ifdef CHECK_RNS
		integer* Acopy=new integer[m*n];
		for(size_t i=0;i<m;i++)
			for(size_t j=0;j<n;j++)
				Acopy[i*n+j]=A[i*lda+j];

#endif
		Givaro::Timer tsplit;tsplit.start();
		size_t mn=m*n;
		double *Arns_tmp = FFLAS::fflas_new<double >(mn*_size*2);
		double* Arns_tmp2= Arns_tmp+mn;
		for (size_t i=0;i<_size;i++)
		  for (size_t j=0;j<mn;j++){
		    uint64_t acci= (uint64_t)Arns[i*rda+j];;
        Arns_tmp [i*2*mn+j] = acci  & ((1<<_shift)-1);
        Arns_tmp2[i*2*mn+j] = (acci >> _shift);
		  }
		
		
#ifdef RNS_DEBUG
		Givaro::ModularExtended<double> ZZ(2UL<<48);;
		std::cout<<"Arns:=";
    FFLAS::WriteMatrix(std::cout, ZZ, _size, mn, Arns, mn);
    //write_field(ZZ,std::cout,Arns, _size, mn,mn,true);
		std::cout<<"Arns1:=";
    FFLAS::WriteMatrix(std::cout, ZZ, _size, mn, Arns_tmp, 2*mn);
    //write_field(ZZ,std::cout,Arns_tmp, _size, mn,2*mn,true);
		std::cout<<"Arns2:=";
    FFLAS::WriteMatrix(std::cout, ZZ, _size, mn, Arns_tmp2, 2*mn);
    //write_field(ZZ,std::cout,Arns_tmp2, _size, mn,2*mn,true);
#endif
		
		integer hM= (_M-1)>>1;
		double *A_beta= FFLAS::fflas_new<double>(2*mn*_ldm);
		double *A_beta2 = A_beta+mn*_ldm;
		tsplit.stop();
#ifdef BENCH_RNS			
		if(m>1 && n>1) std::cerr<<"RNS EXTENDED double (From) - split : "<<tsplit.usertime()<<std::endl;
#endif

		
		Givaro::Timer tfgemmc;tfgemmc.start();
		if (RNS_MAJOR==false)
      // compute A_beta = Ap^T x M_beta
			FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn*2,(int) _ldm,(int) _size, 1.0 , Arns_tmp,(int) mn*2, _crt_out.data(),(int) _ldm, 0., A_beta,(int)_ldm,
                   FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
    //				FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads >());
		
		else // compute A_beta = Ap x M_Beta
			cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)mn*2, (int)_ldm, (int)_size, 1.0 , Arns_tmp, (int)_size, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);

		FFLAS::fflas_delete( Arns_tmp);
		tfgemmc.stop();
#ifdef RNS_DEBUG
		std::cout<<"Abeta1:=";write_field(ZZ,std::cout,A_beta, _ldm, mn,mn,true);
		std::cout<<"Abeta2:=";write_field(ZZ,std::cout,A_beta2, _ldm, mn,mn,true);		
#endif


#ifdef BENCH_RNS			
		if(m>1 && n>1) std::cerr<<"RNS EXTENDED double (From) - fgemm : "<<tfgemmc.usertime()<<std::endl;
#endif
    // compute A using inverse Kronecker transform of A_beta expressed in base 2^log_beta
		integer* Aiter= A;
		size_t k=_ldm;
		size_t k4=((k+3)>>2)+ (((k+3)%4==0)?0:1);
		std::vector<uint16_t> A0(k4<<2,0),A1(k4<<2,0),A2(k4<<2,0),A3(k4<<2,0);
		integer a0,a1,a2,a3,res,res2;
		mpz_t *m0,*m1,*m2,*m3;
		m0= reinterpret_cast<mpz_t*>(&a0);
		m1= reinterpret_cast<mpz_t*>(&a1);
		m2= reinterpret_cast<mpz_t*>(&a2);
		m3= reinterpret_cast<mpz_t*>(&a3);
		mp_limb_t *m0_d,*m1_d,*m2_d,*m3_d;
		m0_d = m0[0]->_mp_d;
		m1_d = m1[0]->_mp_d;
		m2_d = m2[0]->_mp_d;
		m3_d = m3[0]->_mp_d;
		m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = (int) (k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
		m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size  = m3[0]->_mp_size  = (int) (k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
		Givaro::Timer tkroc;
		tkroc.start();
    //		auto sp=SPLITTER();
    //		PARFOR1D(i,m,sp,
		for(size_t i=0;i<m;i++)
			for (size_t j=0;j<n;j++){
				size_t idx=i*n+j;
				for (size_t l=0;l<k;l++){
					uint64_t tmp=(uint64_t)A_beta[l+idx*k];
					uint16_t* tptr= reinterpret_cast<uint16_t*>(&tmp);
					A0[l  ]= tptr[0];
					A1[l+1]= tptr[1];
					A2[l+2]= tptr[2];
					A3[l+3]= tptr[3];
				}
        // see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
				m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A0[0]);
				m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A1[0]);
				m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A2[0]);
				m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A3[0]);
				res = a0;res+= a1;res+= a2;res+= a3;

		
				for (size_t l=0;l<k;l++){
				  uint64_t tmp=(uint64_t)A_beta2[l+idx*k];
					uint16_t* tptr= reinterpret_cast<uint16_t*>(&tmp);
					A0[l  ]= tptr[0];
					A1[l+1]= tptr[1];
					A2[l+2]= tptr[2];
					A3[l+3]= tptr[3];
				}
        // see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
				m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A0[0]);
				m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A1[0]);
				m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A2[0]);
				m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A3[0]);
				res2 = a0;res2+= a1;res2+= a2;res2+= a3;
#ifdef RNS_DEBUG
				std::cout<<"res1:="<<res<<";\n";
				std::cout<<"res2:="<<res2<<";\n";
#endif
          res+=(res2<<_shift);				
				res%=_M;
				
				// get the correct result according to the expected sign of A
				if (res>hM)
					res-=_M;
				if (gamma==0)
					Aiter[j+i*lda]=res;
				else
					if (gamma==integer(1))
						Aiter[j+i*lda]+=res;
					else
						if (gamma==integer(-1))
							Aiter[j+i*lda]=res-Aiter[j+i*lda];
						else{
							Aiter[j+i*lda]*=gamma;
							Aiter[j+i*lda]+=res;
						}

			}
    tkroc.stop();
#ifdef BENCH_RNS			
    if(m>1 && n>1) std::cerr<<"RNS EXTENDED double (From) -  Convert : "<<tkroc.usertime()<<std::endl;
#endif
		m0[0]->_mp_d = m0_d;
		m1[0]->_mp_d = m1_d;
		m2[0]->_mp_d = m2_d;
		m3[0]->_mp_d = m3_d;
		m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc= m3[0]->_mp_alloc = 1;
		m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size = m3[0]->_mp_size  = 0;
		FFLAS::fflas_delete( A_beta);


#ifdef CHECK_RNS
		bool ok=true;
		for (size_t i=0;i<m;i++)
		  for(size_t j=0;j<n;j++)
		    for(size_t k=0;k<_size;k++){
		      int64_t _p =(int64_t) _basis[k];
		      integer curr=A[i*lda+j] - gamma*Acopy[i*n+j];
		      ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);
		      if (!ok) std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
				}
		std::cout<<"RNS EXTENDED double (From) ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif

	}

  

}// end of namespace FFPACK


#endif // __FFLASFFPACK_field_rns_double_extended_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
