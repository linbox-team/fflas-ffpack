/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_pfgemv.inl
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



namespace FFLAS
{
	
	
	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element_ptr
	pfgemv(const Field& F,
		   const FFLAS_TRANSPOSE ta,
		   const size_t m,
		   const size_t n,
		   const typename Field::Element alpha,
		   const typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::ConstElement_ptr X, const size_t incX,
		   const typename Field::Element beta,
		   typename Field::Element_ptr Y, const size_t incY, 
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > & H){
		
		if (H.parseq.numthreads()==1 || m <= 1){
			fgemv(F, ta,  m, n,  alpha, A, lda, X, incX, beta, Y, incY);
			
		}else{
			typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > MMH_t;
			MMH_t H1(H);
			MMH_t H2(H);
			size_t M2 = m>>1;
			PAR_BLOCK{
				
				H1.parseq.set_numthreads(H.parseq.numthreads() >> 1);
				H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());
				
				typename Field::ConstElement_ptr A1 = A;
				typename Field::ConstElement_ptr A2 = A + M2*lda;
				typename Field::Element_ptr C1 = Y;
				typename Field::Element_ptr C2 = Y + M2;
				
				TASK(CONSTREFERENCE(F,H1) MODE( READ(A1,X) READWRITE(C1)),
					 {pfgemv( F, ta,  M2, n, alpha, A1, lda, X, incX, beta, C1, incY, H1);}
					 );
				
				TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,X) READWRITE(C2)),
					 {pfgemv(F, ta, m-M2, n, alpha, A2, lda, X, incX, beta, C2, incY, H2);}
					 );
			}
			
		}
		return Y;		
		
	}
	
	

	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element_ptr
	pfgemv(const Field& F,
		   const FFLAS_TRANSPOSE ta,
		   const size_t m,
		   const size_t n,
		   const typename Field::Element alpha,
		   const typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::ConstElement_ptr X, const size_t incX,
		   const typename Field::Element beta,
		   typename Field::Element_ptr Y, const size_t incY,
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Threads> > & H){
		size_t nt = H.parseq.numthreads(); 	
		if(H.parseq.numthreads()<2|| m <= 1){
			fgemv( F, ta, m, n, alpha, A , lda, X, incX, beta, Y, incY);
		}else{			
			
			if(nt<m){
				size_t N1 = m/nt, N2 = m-N1*nt;
				if(N1>0){
					
					for(size_t i=0; i<N1; i++){
						SYNCH_GROUP(
									FORBLOCK1D(iter, nt,	 H.parseq,  
											   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
													{
														fgemv( F, ta, (iter.end()-iter.begin()), n, alpha, A + iter.begin()*lda, lda, X, incX, beta, Y + iter.begin()*incY, incY);
													} 
													)
											   );
									);
					}
				}//else{}
				if(N2>0){
					SYNCH_GROUP(
								FORBLOCK1D(iter, N2,	 H.parseq,  
										   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
												{
													fgemv( F, ta, (iter.end()-iter.begin()), n, alpha, A + N1*nt + iter.begin()*lda, lda, X, incX, beta, Y + N1*nt + iter.begin()*incY, incY);
												} 
												)
										   );
								);
				}//else{}
				
			}else{
				
				SYNCH_GROUP(
							FORBLOCK1D(iter, m,	 H.parseq,  
									   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
											{
												fgemv( F, ta, (iter.end()-iter.begin()), n, alpha, A + iter.begin()*lda, lda, X, incX, beta, Y + iter.begin()*incY, incY);
											} 
											)
									   );
							);
				
			}
			
		}
		
		return Y;		
		
	}
	
	
	template<class Field, class AlgoT, class FieldTrait>
	typename Field::Element_ptr
	pfgemv(const Field& F,
		   const FFLAS_TRANSPOSE ta,
		   const size_t m,
		   const size_t n,
		   const typename Field::Element alpha,
		   const typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::ConstElement_ptr X, const size_t incX,
		   const typename Field::Element beta,
		   typename Field::Element_ptr Y, const size_t incY,
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > & H){
		
		const int RG = m/H.parseq.numthreads();
		size_t N1 = (m-m/H.parseq.numthreads()*H.parseq.numthreads()>0)?(m-m/H.parseq.numthreads()*H.parseq.numthreads()):0; 
		size_t N2 = H.parseq.numthreads()-N1;
		
		if(H.parseq.numthreads()==1||m<=1){
			fgemv( F, ta, m, n, alpha, A , lda, X, incX, beta, Y, incY);
		}else{
			ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> pH;
			if(N1>0){
				SYNCH_GROUP( 
							
							FORBLOCK1D(iter, N1, pH,   
									   TASK(CONSTREFERENCE(F) MODE( READ(A,X) READWRITE(Y)),
											{
												fgemv( F, ta, RG+1, n, alpha,  A+(RG+1)*iter.blockindex()*lda, lda, X, incX, beta,  Y+iter.blockindex()*(RG+1)*incY, incY);					
											} 
											)
									   );		
							
							 ); 
			}
			
			if(N2>0){
				SYNCH_GROUP( 
							
							FORBLOCK1D(iter, N2, pH,   
									   TASK(CONSTREFERENCE(F) MODE( READ(A,X) READWRITE(Y)),
											{
												
												fgemv( F, ta, RG, n, alpha, A+(RG+1)*N1*lda+RG*iter.blockindex()*lda, lda, X, incX, beta, Y+(RG+1)*N1+iter.blockindex()*RG*incY, incY);		
												
											} 
											)
									   );		
							
							 );  
			}
			
			
		}	
		
		
		return Y;		
		
	}	
	
	template<class Field, class Cut, class Param>
	typename Field::Element_ptr
	pfgemv(const Field& F,
		   const FFLAS_TRANSPOSE ta,
		   const size_t m,
		   const size_t n,
		   const typename Field::Element alpha,
		   const typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::ConstElement_ptr X, const size_t incX,
		   const typename Field::Element beta,
		   typename Field::Element_ptr Y, const size_t incY, 
		   ParSeqHelper::Parallel<Cut,Param> par ){
		
		MMHelper<Field, MMHelperAlgo::Auto, typename FFLAS::ModeTraits<Field>::value, ParSeqHelper::Parallel<Cut,Param> > H (F,m,n,1,par);
		return pfgemv(F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY, H);
	}
	
} // FFLAS

