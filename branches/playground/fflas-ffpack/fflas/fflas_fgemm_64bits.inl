/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
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


uint64_t logtwoceil(uint64_t v) {
	uint64_t r(0);
	while (v >>= 1) ++r;
	return ++r;
}

void fgemm_mp_kron2(int64_t p, size_t m, size_t n, size_t k, 
		    int64_t alpha, const int64_t* A, size_t lda, const int64_t* B, size_t ldb,
		    int64_t beta, int64_t* C, size_t ldc);

void fgemm_mp_kron3(int64_t p, size_t m, size_t n, size_t k, 
		    int64_t alpha, const int64_t* A, size_t lda, const int64_t* B, size_t ldb,
		    int64_t beta, int64_t* C, size_t ldc);

// C= alpha.AB+beta.C
void fgemm(int64_t p, size_t m, size_t n, size_t k, 
	   int64_t alpha,  const int64_t* A, size_t lda, const int64_t* B, size_t ldb,
	   int64_t beta, int64_t* C, size_t ldc)
{
	uint64_t logp,logk;
	logp=logtwoceil(p); 
	logk=logtwoceil(k);
	
	if (logp+logk<53){ // use Karatsuba
		fgemm_mp_kron2(p,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
		return;
	}
	if (2*logp+3*logk<159){ // use Toom-3
		fgemm_mp_kron3(p,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
		return;
	}
}




struct Kronecker_modp {
	uint64_t           p;
	uint64_t        logp;
	size_t      log_beta;
	uint64_t        mask;
	double            nu;
	

	Kronecker_modp(uint64_t prime, size_t k) 
		: p(prime), logp(logtwoceil(p)), log_beta(logp/k),
		  mask( (1UL<<log_beta)-1UL)
	{
		nu=(1UL<<logp)*((double)(1UL<<log_beta)/(double)p);
		nu/=(double)(1UL<<log_beta);
	}
	
	inline void split(int64_t  a, double& u){u= a;}
	
	template<typename... T>
	inline void split(int64_t  a, double& x, T&... u){
		x= a&mask;
		split(a>>log_beta,u...);
	}
	
	inline void reconstruct(uint64_t& r, double x){r = (uint64_t)x %p;}
	
	template<typename... T>
	inline void reconstruct(uint64_t& r, double x, T... u){
		reconstruct(r,u...);
		r<<= log_beta;
		r+=(uint64_t)x;
		r%=p;
	}

	inline void reconstruct_over(uint64_t& r, double x){r = (uint64_t)x %p;}

	template<typename... T>
	inline void reconstruct_over(uint64_t& r, double x, T... u){
		reconstruct_over(r,u...);
		uint64_t q;
		double c1;
		c1= (double)(r>>(logp-log_beta));
		q= (uint64_t)(c1*nu);
		r&=((1UL<<(64-log_beta))-1);
		r= (r<<log_beta) +(uint64_t)x;
		r-=q*p; 
		while (r>p) r-=p;
	}

	bool overflow(size_t n)const {return (log_beta+logp) >64;}

};

void fgemm_mp_kron2(int64_t p, size_t m, size_t n, size_t k, 
		    int64_t alpha, const int64_t* A, size_t lda,const int64_t* B, size_t ldb,
		    int64_t beta, int64_t* C, size_t ldc){
	double *Ad0, *Bd0, *Cd0, *Ad1, *Bd1, *Cd1, *Cd01;
	size_t mn,mk,kn;
	if (alpha!=1 && alpha!=(p-1) && alpha!= -1){ 
		cout<<"alpha != 1 or -1 in fgemm_mp_kron2 (not yet supported)...aborting"<<endl;
		cout<<"alpha="<<alpha<<endl;
		cout<<"p="<<p<<endl;
		exit(1);
	}
		
	mn=m*n;
	mk=m*k;
	kn=k*n;
	Ad0  = new double[mk];
	Ad1  = new double[mk];
	Bd0  = new double[kn];
	Bd1  = new double[kn];
	Cd0  = new double[mn];
	Cd1  = new double[mn];
	Cd01 = new double[mn];

	Kronecker_modp Kron(p[0],2);
	for(size_t i = 0,t=0; i <m; ++i) 
		for(size_t j=0;j<k;j++,t++)
			Kron.split(A[i*ldA+j],Ad0[t],Ad1[t]);
	for(size_t i = 0,t=0; i <k; ++i) 
		for(size_t j=0;j<n;j++,t++)
			Kron.split(B[i*ldb+j],Bd0[t],Bd1[t]);

	FFLAS::DoubleDomain DD;    
	double aaa=(alpha==1)?1:1;
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, aaa, Ad0, k, Bd0, n, DD.zero, Cd0, n);
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, k, aaa, Ad1, k, Bd1, n, DD.zero, Cd1, n);

	for(size_t i = 0; i <mk; ++i) 
		Ad0[i] -= Ad1[i];
	for(size_t i = 0; i <kn; ++i) 
		Bd1[i] -= Bd0[i];
    
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, aaa, Ad0, n, Bd1, n, DD.zero, Cd01, n);

	for(size_t i = 0; i <mn; ++i) {
		Cd01[i] += Cd0[i];
		Cd01[i] += Cd1[i];
	}

	uint64_t tmp;
	if (Kron.overflow(k)){//cout<<"USING PBR"<<endl;
		for(size_t i = 0,t=0; i <m; ++i) 
			for(size_t j=0;j<n;j++,t++){
				Kron.reconstruct_over(tmp, Cd0[t], Cd01[t], Cd1[t]);
				if (beta==0) C[i*ldc+j]=tmp;				
				else
					if (beta==1) C[i*ldc+j]+=tmp;
					else if (beta==-1) C[i*ldc+j]-=tmp;
					else C[i*ldc+j]=tmp+beta*C[i*ldc+j];
			}
	}
	else{
		for(size_t i = 0,t=0; i <m; ++i) 
			for(size_t j=0;j<n;j++,t++){
				Kron.reconstruct(tmp, Cd0[t], Cd01[t], Cd1[t]);
				if (beta==0) C[i*ldc+j]=tmp;				
				else
					if (beta==1) C[i*ldc+j]+=tmp;
					else if (beta==-1) C[i*ldc+j]-=tmp;
					else C[i*ldc+j]=tmp+beta*C[i*ldc+j];
			}
	}
	
	delete [] Ad0;
	delete [] Bd0;
	delete [] Cd0;
	delete [] Ad1;
	delete [] Bd1;
	delete [] Cd1;
	delete [] Cd01;
}

void fgemm_mp_kron3(int64_t p, size_t m, size_t n, size_t k,
		    int64_t alpha, const int64_t* A, size_T lda, const int64_t* B, size_t ldb,
		    int64_t beta, int64_t* C, size_t ldc){
	double *Ad0, *Ad1, *Ad2, *Ad02, *Bd0, *Bd1, *Bd2, *Bd02, *Cd0, *Cd1, *Cd2, *Cd3, *Cd4;
	size_t mk=m*k,kn=k*n,mn=m*n;	
	if (alpha!=1 ){
		cout<<"alpha != 1 in fgemm_mp_kron3 (not yet supported)...aborting"<<endl;
		cout<<"alpha="<<alpha<<endl;
		cout<<"p="<<p<<endl;
		exit(1);
	}

	Ad0  = new double[mk];
	Ad1  = new double[mk];
	Ad2  = new double[mk];
	Ad02 = new double[mk];
	Bd0  = new double[kn];
	Bd1  = new double[kn];
	Bd2  = new double[kn];
	Bd02 = new double[kn];
	Cd0  = new double[mn];
	Cd1  = new double[mn];
	Cd2  = new double[mn];
	Cd3  = new double[mn];
	Cd4  = new double[mn];
	double T1,T2;
    
	Kronecker_modp Kron(p[0],3);
	for(size_t i = 0,t=0; i <m; ++i) 
		for(size_t j=0;j<k;j++,t++)
			Kron.split(A[i*ldA+j],Ad0[t],Ad1[t],Ad2[t]);
	for(size_t i = 0,t=0; i <k; ++i) 
		for(size_t j=0;j<n;j++,t++)
			Kron.split(B[i*ldb+j],Bd0[t],Bd1[t],Bd2[t]);

	
	FFLAS::DoubleDomain DD;    
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, DD.one, Ad0, k, Bd0, n, DD.zero, Cd0, n);
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, DD.one, Ad2, k, Bd2, n, DD.zero, Cd4, n);

	for(size_t i = 0; i <mk; ++i) {
		Ad02[i] = Ad0[i]+Ad2[i];
		Ad0[i] += (Ad1[i]*2.)+(Ad2[i]*4.);
	}
	for(size_t i = 0; i <kn; ++i) {
		Bd02[i] = Bd0[i]+Bd2[i];
		Bd0[i] += (Bd1[i]*2.)+(Bd2[i]*4.);
	}
    
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, DD.one, Ad0, k, Bd0, n, DD.zero, Cd3, n);

	for(size_t i = 0; i <mk; ++i) {
		Ad0[i] =Ad02[i]-Ad1[i];
		Ad1[i]+=Ad02[i];
	}
	for(size_t i = 0; i <kn; ++i) {
		Bd0[i] =Bd02[i]-Bd1[i];
		Bd1[i]+=Bd02[i];		
	}

	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, DD.one, Ad0, k, Bd0, n, DD.zero, Cd2, n);
	FFLAS::fgemm (DD, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, k, DD.one, Ad1, k, Bd1, n, DD.zero, Cd1, n);

	for(size_t i = 0; i <mn; ++i) {
		T1=(3.*Cd0[i]+(2.*Cd2[i])+Cd3[i])/6. - 2.*Cd4[i] ;
		T2=(Cd1[i]+Cd2[i])/2.;
		Cd1[i]-=T1;
		Cd2[i]=T2-Cd0[i]-Cd4[i];
		Cd3[i]=T1-T2;
	}



	uint64_t tmp;
	if (Kron.overflow(k)){cout<<"USING PBR"<<endl;
		for(size_t i = 0,t=0; i <m; ++i) 
			for(size_t j=0;j<n;j++,t++){
				Kron.reconstruct_over(tmp, Cd0[t], Cd1[t], Cd2[t], Cd3[t], Cd4[t]);
				if (beta==0) C[i*ldc+j]=tmp;				
				else
					if (beta==1) C[i*ldc+j]+=tmp;
					else if (beta==-1) C[i*ldc+j]-=tmp;
					else C[i*ldc+j]=tmp+beta*C[i*ldc+j];
			}
	}
	else{
		for(size_t i = 0,t=0; i <m; ++i) 
			for(size_t j=0;j<n;j++,t++){
				Kron.reconstruct(tmp, Cd0[i], Cd1[i], Cd2[i], Cd3[i], Cd4[i]);
				if (beta==0) C[i*ldc+j]=tmp;				
				else
					if (beta==1) C[i*ldc+j]+=tmp;
					else if (beta==-1) C[i*ldc+j]-=tmp;
					else C[i*ldc+j]=tmp+beta*C[i*ldc+j];
			}
	}
	delete [] Ad0;
	delete [] Ad1;
	delete [] Ad2;
	delete [] Ad02;
	delete [] Bd0;
	delete [] Bd1;		
	delete [] Bd2;
	delete [] Bd02;	
	delete [] Cd0;
	delete [] Cd1;	
	delete [] Cd2;
	delete [] Cd3;
	delete [] Cd4;
	

}

