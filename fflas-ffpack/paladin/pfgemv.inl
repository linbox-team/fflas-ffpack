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


		if (m<2){
/*
SYNCH_GROUP(
FFLAS::WriteMatrix (std::cout << "A:in ="<< std::endl, F, m, n, A, lda) << std::endl;
FFLAS::WriteMatrix (std::cout << "X:in ="<< std::endl, F, n, incX, X, incX) << std::endl;
WAIT;
);
*/
			/*
			  PAR_BLOCK{
			  SYNCH_GROUP(
							TASK(CONSTREFERENCE(F) MODE( READ(A,X) READWRITE(Y)),
							{fgemv(F, ta,  m, n,  alpha, A, lda, X, incX, beta, Y, incY);} 
							);
							);
							}
			*/
			fgemv(F, ta,  m, n,  alpha, A, lda, X, incX, beta, Y, incY);
	 
	//FFLAS::WriteMatrix (std::cout << "Y:out ="<< std::endl, F, incY, incY, Y, incX) << std::endl;
		}else{
			typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > MMH_t;
			MMH_t H1(H);
			MMH_t H2(H);
			size_t M2 = m>>1;
			PAR_BLOCK{
SYNCH_GROUP(
				if(H1.parseq.numthreads()>1){
					H1.parseq.set_numthreads(H.parseq.numthreads() >> 1);
					H2.parseq.set_numthreads(H.parseq.numthreads() - H1.parseq.numthreads());
				}else{
					H1.parseq.set_numthreads(H.parseq.numthreads());
					H2.parseq.set_numthreads(H.parseq.numthreads());	
				}
				
				typename Field::ConstElement_ptr A1 = A;
				typename Field::ConstElement_ptr A2 = A + M2*lda;
				typename Field::Element_ptr C1 = Y;
				typename Field::Element_ptr C2 = Y + M2;
//FFLAS::WriteMatrix (std::cout << "A1:in ="<< std::endl, F, M2, n, A1, incX) << std::endl;	
				// 2 multiply (1 split on dimension m)	
	//std::cout<<"m1: "<<M2<<std::endl;
				TASK(CONSTREFERENCE(F,H1) MODE( READ(A1,X) READWRITE(C1)),
					 {pfgemv( F, ta,  M2, n, alpha, A1, lda, X, incX, beta, C1, incY, H1);}
					 );
//WAIT;
	//std::cout<<"m2: "<<m-M2<<std::endl;
//FFLAS::WriteMatrix (std::cout << "A2:in ="<< std::endl, F, m-M2, n, A2, incX) << std::endl;	
				TASK(MODE(CONSTREFERENCE(F,H2) READ(A2,X) READWRITE(C2)),
					 {pfgemv(F, ta, m-M2, n, alpha, A2, lda, X, incX, beta, C2, incY, H2);}
					 );
			}
//WAIT;
);
		}
		return Y;		
		
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		
		typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Threads> > MMH_t;
		MMH_t HH(H);
		size_t BS;
		if(BS<H.parseq.numthreads()||BS==0) BS = m/(H.parseq.numthreads());
		size_t Niter = std::max(m, H.parseq.numthreads());
//std::cout<<" H.parseq.numthreads(): "<< H.parseq.numthreads()<<std::endl;
		if(m<Niter){

SYNCH_GROUP(

			PAR_BLOCK{
				
				HH.parseq.set_numthreads(H.parseq.numthreads());
				using FFLAS::CuttingStrategy::Row;
				using FFLAS::StrategyParameter::Threads;				
				
				FORBLOCK1D(iter, m,	SPLITTER(BS,Row,Threads),  
						   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(YY)),
								{
fgemv( F, ta, 1, n, alpha, A + iter.blockindex()*lda, lda, X, incX, beta, Y+ iter.blockindex(), incY);
								} 
							)
				);
			}
);
			
		}else{


				if(H.parseq.numthreads()>1){
			PAR_BLOCK{

				HH.parseq.set_numthreads(H.parseq.numthreads());
				using FFLAS::CuttingStrategy::Row;
				using FFLAS::StrategyParameter::Threads;				
				
				FORBLOCK1D(iter, m,	SPLITTER(BS,Row,Threads), 		
						   TASK(CONSTREFERENCE(F,H) MODE( READ(A1,X) READWRITE(YY)),
								{
//FFLAS::WriteMatrix (std::cout << "A1:in ="<< std::endl, F, 1, n, A + iter.blockindex()*(H.parseq.numthreads()-iter.blockindex())*lda, lda) << std::endl;	
fgemv( F, ta, 1, n, alpha, A + iter.blockindex()*(H.parseq.numthreads()-iter.blockindex())*lda, lda, X, incX, beta, Y+ iter.blockindex()*(H.parseq.numthreads()-iter.blockindex())*incY, incY);
}
								)
						   
						   );
			}
				}else{
			PAR_BLOCK{

				HH.parseq.set_numthreads(H.parseq.numthreads());
				using FFLAS::CuttingStrategy::Row;
				using FFLAS::StrategyParameter::Threads;				
				
				FORBLOCK1D(iter, m,	SPLITTER(BS,Row,Threads), 		
						   TASK(CONSTREFERENCE(F,H) MODE( READ(A1,X) READWRITE(YY)),
								{
//FFLAS::WriteMatrix (std::cout << "A1:in ="<< std::endl, F, 1, n, A + iter.blockindex()*lda, lda) << std::endl;WAIT;
fgemv( F, ta, 1, n, alpha, A + iter.blockindex()*lda, lda, X, incX, beta, Y + + iter.blockindex()*incY, incY);
}
								)
						   
						   );
			}
				}

		}


		return Y;		
		
	}
	
/*
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
		   size_t GS_Cache,
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > & H){
		
		typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > MMH_t;
		MMH_t HH(H);

		//if m/GS_Cache>0 then split each line of the matrix into m/GS_Cache blocks
		size_t NbBlocs;
		if(GS_Cache>0) NbBlocs = n/GS_Cache; 
		//else split the matrix into m/(H.parseq.numthreads()) line-wise blocks
		else  NbBlocs = n/(H.parseq.numthreads());
		//distribute each block into one thread
		if(NbBlocs>0){

			PAR_BLOCK{
				for(size_t j=0;j<m;j++) 
					{

						HH.parseq.set_numthreads(H.parseq.numthreads());
						using FFLAS::CuttingStrategy::Row;
						using FFLAS::StrategyParameter::Grain;		
						
						//for(size_t iter=0;iter<NbBlocs;iter++)
						FORBLOCK1D(iter, NbBlocs,	SPLITTER(GS_Cache,Row,Grain), 			
								   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
										{
											fgemv( F, ta, 1, GS_Cache, alpha, (A+j) + iter.blockindex()*GS_Cache, GS_Cache, X+ iter.blockindex()*GS_Cache, incX, beta, Y+j, incY);
											fadd (F, 1, GS_Cache, Y+j, incY, Y+j, incY, Y+j, incY);
										}
										)
								   );
						if(n - NbBlocs*n>0){
								   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
										{
											fgemv( F, ta, 1, GS_Cache, alpha, (A+j) + (n - NbBlocs), (n - NbBlocs), X+ (n - NbBlocs), incX, beta, Y+j, incY);
											fadd (F, 1, GS_Cache, Y+j, incY, Y+j, incY, Y+j, incY);
										}

								   );
						}
					}
			}

		}else{
			size_t Niter = std::max(m, H.parseq.numthreads());
			if(m<Niter){
				PAR_BLOCK{
				
					HH.parseq.set_numthreads(H.parseq.numthreads());
					using FFLAS::CuttingStrategy::Row;
					using FFLAS::StrategyParameter::Grain;				
					
					FORBLOCK1D(iter, Niter,	SPLITTER(GS_Cache,Row,Grain),			
							   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
									{fgemv( F, ta, 1, n, alpha, A + iter.blockindex()*lda, lda, X, incX, beta, Y+ iter.blockindex(), incY);}
									)
							   );
				}
						
			}else{
				PAR_BLOCK{
				
					HH.parseq.set_numthreads(H.parseq.numthreads());
					using FFLAS::CuttingStrategy::Row;
					using FFLAS::StrategyParameter::Grain;				
					
					FORBLOCK1D(iter, Niter,	SPLITTER(GS_Cache,Row,Grain), 			
							   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
									{fgemv( F, ta, 1, n, alpha, A + iter.blockindex()*GS_Cache*lda, lda, X, incX, beta, Y+ iter.blockindex(), incY);}
									)
							   );
				}
			}
		}
		
		
		return Y;		
		
	}
*/
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
		   size_t GS_Cache,
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > & H){
		
		typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > MMH_t;
		MMH_t HH(H);

		size_t N = min(m,n);
		const int TILE = min(min(m,GS_Cache), min(n,GS_Cache) ); 
		//Compute tiles in each dimension
		const int nEven = N - N%TILE;
		//const int wTiles = nEven/TILE;

		PAR_BLOCK{

			//Main body of the matrix
			/*
			Potential loop coalescing to improve data locality
			Then possible to put into FORBLOCK1D(divide the taks as even as possible)
			*/
			for(int ii=TILE; ii<nEven; ii+=TILE){
				for(int jj=TILE; jj<nEven; jj+=TILE){
					fgemv( F, ta, TILE, TILE, alpha, A+ii*TILE+jj, TILE, X+ii*TILE, incX, beta, Y+ii*TILE, incY);
				}
			}

			if(N-nEven>0){
					//The bottom rows of the matrix
					fgemv( F, ta, N-nEven, n, alpha, A+(N-nEven), N-nEven, X, incX, beta, Y+(N-nEven), incY);

					//The right columns of the matrix
					fgemv( F, ta, m, N-nEven, alpha, A, m, X, incX, beta, Y, incY);

					}
		
					//The bottom right corner of the matrix
					fgemv( F, ta, N-nEven, N-nEven, alpha, A+nEven*lda+(N-nEven), N-nEven, X+N-nEven, incX, beta, Y+N-nEven, incY);
		}
		return Y;		
		
	}	
	
} // FFLAS

