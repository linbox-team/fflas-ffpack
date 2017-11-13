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

#include  "mpi.h"

#include "fflas-ffpack/paladin/pfgemv.inl"




template <typename T > class chooseMPItype;
template <> struct chooseMPItype<double>{ static constexpr MPI_Datatype val = MPI_DOUBLE;};
template <> struct chooseMPItype<float>{ static constexpr MPI_Datatype val = MPI_FLOAT;};
template <> struct chooseMPItype<int32_t>{ static constexpr MPI_Datatype val = MPI_INT32_T;};
template <> struct chooseMPItype<int64_t>{ static constexpr MPI_Datatype val = MPI_INT64_T;};






void printMPItype(MPI_Datatype type){ 
	if(MPI_DOUBLE==type)  std::cout<<"MPI_DOUBLE is used"<< std::endl;
	if(MPI_FLOAT==type)   std::cout<<"MPI_FLOAT is used"<< std::endl;
	if(MPI_INT32_T==type) std::cout<<"MPI_INT32_T is used"<<std::endl;
	if(MPI_INT64_T==type) std::cout<<"MPI_INT64_T is used"<<std::endl;
};






namespace FFLAS
{
	
	
	
	template<typename  Field, class AlgoT, class FieldTrait>
	typename Field::Element_ptr
	pfgemv_mpi(const Field& F,
			   const FFLAS_TRANSPOSE ta,
			   const size_t m,
			   const size_t n,
			   const typename Field::Element alpha,
			   const typename Field::ConstElement_ptr A, const size_t lda,
			   const typename Field::ConstElement_ptr X, const size_t incX,
			   const typename Field::Element beta,
			   typename Field::Element_ptr Y, const size_t incY, 
			   MMHelper<Field, AlgoT, FieldTrait, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> > & H){
		
		int	rank, nprocs;
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Status status[nprocs];
		MPI_Request recv_request[nprocs];
		
		if(rank==0) 
			{	std::cout << "<<<<<<<<<<<Master thread: ";printMPItype(chooseMPItype<typename Field::Element>::val);	
				FFLAS::WriteMatrix (std::cout << "A:="<<std::endl, F, m, m, A, lda) << std::endl;
				FFLAS::WriteMatrix (std::cout << "X:="<<std::endl, F, m, incX, X, incX) << std::endl;
				FFLAS::WriteMatrix (std::cout << "Y:="<<std::endl, F, m, incY, Y, incY) << std::endl;
				
				if(nprocs==1) {
					PAR_BLOCK {
						pfgemv(F,ta,m, n, alpha, A, lda, X,incX,beta,Y, incY, H);
					}
					return Y;
				}else{

				std::cerr<<"Master: typename Field::Element = "<<typeid(typename Field::Element).name()<<std::endl;
	if(m%nprocs!=0){
						
						for(int i=0;i<m%nprocs;i++){
							
							MPI_Isend(A+i*(m/nprocs+1)*lda, (m/nprocs+1)*n,  chooseMPItype<typename Field::Element>::val, i+1, 123, MPI_COMM_WORLD, &recv_request[i+1]);
							//FFLAS::WriteMatrix (std::cout << "A sent to proc("<<i+1<<"):"<<std::endl, F, (m/nprocs+1), n, A+i*(m/nprocs+1)*lda, lda) << std::endl;
						}
						for(int i=0;i<nprocs-m%nprocs-1;i++){
							
							MPI_Isend(A+(m%nprocs)*(m/nprocs+1)*lda+i*(m/nprocs)*lda, (m/nprocs)*n,  chooseMPItype<typename Field::Element>::val, i+m%nprocs+1, 123, MPI_COMM_WORLD, &recv_request[i+m%nprocs+1]);
							//FFLAS::WriteMatrix (std::cout << "A sent to proc("<<i+m%nprocs+1<<"):"<<std::endl, F, (m/nprocs), n, A+(m%nprocs)*(m/nprocs+1)*lda+i*(m/nprocs)*lda, lda) << std::endl;
						}
						
						
						for(int i=1;i<nprocs;i++){
							
							MPI_Isend(X, n,  chooseMPItype<typename Field::Element>::val, i, 123, MPI_COMM_WORLD, &recv_request[i]);
							//FFLAS::WriteMatrix (std::cout << "X sent to proc("<<i+1<<"):"<<std::endl, F, n, incX, X, incX) << std::endl;
						}
						
						std::cerr<<"Master after send  : typename Field::Element = "<<typeid(typename Field::Element).name()<<std::endl;

						
						
						
						PAR_BLOCK {
							pfgemv(F,ta,m-(m%nprocs)*(m/nprocs+1)-(nprocs-m%nprocs-1)*(m/nprocs), n, alpha, A+(m%nprocs)*(m/nprocs+1)*lda+(nprocs-m%nprocs-1)*(m/nprocs)*lda, lda, X,incX,beta,Y+(m%nprocs)*(m/nprocs+1)*incY+(nprocs-m%nprocs-1)*(m/nprocs)*incY, incY, H);
						}

						
						for(int i=0;i<m%nprocs;i++){
							MPI_Recv(Y+i*(m/nprocs+1)*incY, (m/nprocs+1)*incY,  chooseMPItype<typename Field::Element>::val, i+1, 124, MPI_COMM_WORLD, &status[i+1]);
							//FFLAS::WriteMatrix (std::cout << "YY received from proc("<<i+1<<"):"<<std::endl, F, (m/nprocs+1)*incY, incY, Y, incY) << std::endl;
						}
						
						for(int i=0;i<nprocs-m%nprocs-1;i++){
							MPI_Recv(Y+(m%nprocs)*(m/nprocs+1)*incY+i*(m/nprocs)*incY, (m/nprocs)*incY,  chooseMPItype<typename Field::Element>::val, i+m%nprocs+1, 124, MPI_COMM_WORLD, &status[i+m%nprocs+1]);
							//FFLAS::WriteMatrix (std::cout << "YY received from proc("<<i+m%nprocs+1<<"):"<<std::endl, F, (m/nprocs)*incY, incY, Y, incY) << std::endl;
						}
						
						
						
						
						
					} else{ //m%nprocs==0
						
						for(int i=0;i<nprocs-1;i++){
							
							MPI_Isend(A+i*(m/nprocs)*lda, (m/nprocs)*n,  chooseMPItype<typename Field::Element>::val, i+1, 123, MPI_COMM_WORLD, &recv_request[i+1]);
							
						}
						
						for(int i=1;i<nprocs;i++){
							
							MPI_Isend(X, n,  chooseMPItype<typename Field::Element>::val, i, 123, MPI_COMM_WORLD, &recv_request[i]);
							
						}
						
						
						
						
						
						
						PAR_BLOCK {
							pfgemv(F,ta,m-(nprocs-1)*(m/nprocs), n, alpha, A+(nprocs-1)*(m/nprocs)*lda, lda, X,incX,beta,Y+(nprocs-1)*(m/nprocs)*incY, incY, H);
						}
						
						for(int i=0;i<nprocs-1;i++){
							
							MPI_Recv(Y+i*(m/nprocs)*incY, (m/nprocs)*incY,  chooseMPItype<typename Field::Element>::val, i+1, 124, MPI_COMM_WORLD, &status[i+1]);
							//FFLAS::WriteMatrix (std::cout << "YY received from proc("<<i+1<<"):"<<std::endl, F, (m/nprocs)*incY, incY, Y, incY) << std::endl;
						}
						
						
						
					}
					
					

					
					
				}

				//tag=false;
				//for(int i=1; i<nprocs; i++)MPI_Isend(&tag, 1, MPI_C_BOOL, i, 123, MPI_COMM_WORLD, &recv_request[i]);				

				return Y;
				
				
			} else { //rank!=0
				
			
			//while(1){
				

				typename Field::Element_ptr AA, XX, YY;
				//if(tag){
std::cout << ">>>>>>>>>>Worker thread: ";printMPItype(chooseMPItype<typename Field::Element>::val);
std::cerr<<"typename Field::Element = "<<typeid(typename Field::Element).name()<<std::endl;
					if(rank<m%nprocs+1){
						typename Field::Element_ptr AA = fflas_new(F,m/nprocs+1,n);
						MPI_Recv(AA, (m/nprocs+1)*n,  chooseMPItype<typename Field::Element>::val, 0, 123, MPI_COMM_WORLD,&status[0]);
						
						
						XX = fflas_new(F,n,incX);
						
						MPI_Recv(XX, n,  chooseMPItype<typename Field::Element>::val, 0, 123, MPI_COMM_WORLD,&status[0]);
						
						
						YY = fflas_new(F,m/nprocs+1, incY);
						PAR_BLOCK {
							
							//FFLAS::WriteMatrix (std::cout <<"rank("<<rank<<"): A received from proc(0):"<<std::endl, F, m/nprocs+1, n, AA, lda) << std::endl;
							//FFLAS::WriteMatrix (std::cout <<"rank("<<rank<<"): X received from proc(0):"<<std::endl, F, n, incX, XX, incX) << std::endl;
							pfgemv(F,ta,m/nprocs+1, n, alpha, AA, lda, XX, incX, beta, YY, incY, H);

							
						}
						FFLAS::fflas_delete(AA);
						FFLAS::fflas_delete(XX);				
						
						MPI_Isend(YY, (m/nprocs+1),  chooseMPItype<typename Field::Element>::val, 0, 124, MPI_COMM_WORLD, &recv_request[rank]);
						FFLAS::fflas_delete(YY);
						
					}else{ 
						
						AA = fflas_new(F,m/nprocs,n);
						MPI_Recv(AA, (m/nprocs)*n,  chooseMPItype<typename Field::Element>::val, 0, 123, MPI_COMM_WORLD,&status[0]);
						
						
						XX = fflas_new(F,n,incX);
						
						MPI_Recv(XX, n,  chooseMPItype<typename Field::Element>::val, 0, 123, MPI_COMM_WORLD,&status[0]);
						
						
						YY = fflas_new(F,m/nprocs,incY);
						
						PAR_BLOCK {
							
							//FFLAS::WriteMatrix (std::cout <<"rank("<<rank<< "): A received from proc(0):"<<std::endl, F, m/nprocs, n, AA, lda) << std::endl;
							//FFLAS::WriteMatrix (std::cout <<"rank("<<rank<<"): X received from proc(0):"<<std::endl, F, n, incX, XX, incX) << std::endl;
							pfgemv(F,ta,m/nprocs, n, alpha, AA, lda, XX, incX, beta, YY, incY, H);

						}
						
						FFLAS::fflas_delete(AA);
						FFLAS::fflas_delete(XX);		

						MPI_Isend(YY, (m/nprocs),  chooseMPItype<typename Field::Element>::val, 0, 124, MPI_COMM_WORLD, &recv_request[rank]);
						
						FFLAS::fflas_delete(YY);						

					}
				//}//if(tag)
MPI_Isend(Y, (m/nprocs),  chooseMPItype<typename Field::Element>::val, 0, 124, MPI_COMM_WORLD, &recv_request[rank]);
			//}//while(1)
			
		}   
		

		
	}
	
	
	
	
	
} // FFLAS

