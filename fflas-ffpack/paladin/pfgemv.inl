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
std::cout<<"H.parseq.numthreads():"<<H.parseq.numthreads()<<std::endl;

		if (H.parseq.numthreads()<2){
			fgemv(F, ta,  m, n,  alpha, A, lda, X, incX, beta, Y, incY);

		}else{
			typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > MMH_t;
			MMH_t H1(H);
			MMH_t H2(H);
			size_t M2 = m>>1;
			PAR_BLOCK{

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

		if(H.parseq.numthreads()<2){
	 fgemv( F, ta, m, n, alpha, A , lda, X, incX, beta, Y, incY);
}else{

//SYNCH_GROUP(

			PAR_BLOCK{
				using FFLAS::CuttingStrategy::Row;
				using FFLAS::StrategyParameter::Threads;				
typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Threads> > MMH_t;
MMH_t HH(H);
ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Threads> pH;

				HH.parseq.set_numthreads(H.parseq.numthreads());

				FORBLOCK1D(iter, m,	 pH,  
						   TASK(CONSTREFERENCE(F) MODE( READ(A1,X) READWRITE(Y)),
								{
 pfgemv( F, ta, (iter.end()-iter.begin()), n, alpha, A + iter.begin()*lda, lda, X, incX, beta, Y + iter.begin(), incY, HH);
								} 
							)
				);
			}
//);
			

}

		return Y;		
		
	}
	

/////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Field>
	void partfgemv(const Field& F,
		   const FFLAS_TRANSPOSE ta,
		   const size_t m,
		   const size_t n,
		   const typename Field::Element alpha,
		   const typename Field::ConstElement_ptr A, const size_t lda,
		   const typename Field::ConstElement_ptr X, const size_t incX,
		   const typename Field::Element beta,
		   typename Field::Element_ptr Y, const size_t incY){

		typename Field::Element_ptr tmp  = FFLAS::fflas_new(F, n, incY);
		fassign (F, n, incY, Y, incY, tmp, incY);
		fgemv( F, ta, m, n, alpha, A, lda, X, incX, beta, Y, incY);
		fadd (F, n, incY, tmp, incY, alpha, Y, incY, Y, incY);

		return;
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
		   size_t GS_Cache,
		   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > & H){



FFLAS::WriteMatrix (std::cout << "A:in ="<< std::endl, F, m, n, A, lda) << std::endl;
FFLAS::WriteMatrix (std::cout << "X:in ="<< std::endl, F, m, 1, X, incX) << std::endl;
FFLAS::WriteMatrix (std::cout << "Y:in ="<< std::endl, F, m, 1, Y, incY) << std::endl;



		typedef MMHelper<Field,AlgoT,FieldTrait,ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Grain> > MMH_t;
		MMH_t HH(H);

		size_t N = min(m,n);
		const int TILE = min(min(m,GS_Cache), min(n,GS_Cache) ); 
		//Compute tiles in each dimension
		const int nEven = N - N%TILE;
		const int wTiles = nEven/TILE;





SYNCH_GROUP(		
PAR_BLOCK{
std::cout<<" H.parseq.numthreads(): "<< H.parseq.numthreads()<<std::endl;
			//Main body of the matrix

			for(int ii=0; ii<nEven; ii+=TILE){
				for(int jj=0; jj<ii; jj+=TILE){
//-----------------------------------------------------------------------------------------------------------
			for(int j=jj; j<jj+TILE; j++){
				for(int i=ii; i<ii+TILE; i++){
FFLAS::WriteMatrix (std::cout << "A1:in ="<< std::endl, F, TILE, TILE, A+i*lda+j, lda) << std::endl;	
FFLAS::WriteMatrix (std::cout << "X1:in ="<< std::endl, F, TILE, 1, X+j, incX) << std::endl;
FFLAS::WriteMatrix (std::cout << "Y1:in ="<< std::endl, F, TILE, 1, Y+j, incY) << std::endl;
					//fgemv( F, ta, TILE, TILE, alpha, A+ii*lda+jj, lda, X+ii*TILE, incX, beta, Y+jj, incY);
					//partfgemv( F, ta, TILE, TILE, alpha, A+ii*lda+jj, lda, X+ii*TILE, incX, beta, Y+jj, incY);
FFLAS::WriteMatrix (std::cout << "Y1:out ="<< std::endl, F, TILE, 1, Y+j, incY) << std::endl;
FFLAS::WriteMatrix (std::cout << "Y2:out ="<< std::endl, F, TILE, 1, Y+j, incY) << std::endl;
				}
			}  
//-----------------------------------------------------------------------------------------------------------
				}
			}  

		}   ); //SYNCH_GROUP

		return Y;		
		
	}	
	
} // FFLAS

