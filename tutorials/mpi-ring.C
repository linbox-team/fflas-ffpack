#include  "mpi.h"
#include <unistd.h>
#include "stdlib.h"
#include<iostream>

#include <fflas-ffpack/fflas/fflas.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/fflas_io.h"
#include <fflas-ffpack/ffpack/ffpack.h>


#define MAXSIZE 2


using Givaro::Modular;
using namespace FFLAS;
using namespace FFPACK;


template<class Field>
void unpack( const Field &F, typename Field::Element_ptr A,  const size_t M, const size_t N, const size_t lda, double data[]){

 
  for(size_t i=0; i< M; i++){
    for(size_t j=0; j< N; j++){
        data[i*lda+j]=A[i*lda+j];
    }
  }

}

template<typename Field>
void repack(const Field &F, typename Field::Element_ptr A, const size_t M, const size_t N, const size_t lda, double data[]){
  for(size_t i=0; i< M; i++){
    for(size_t j=0; j< N; j++){
        F.assign(*(A+i*lda+j), data[i*lda+j]);
    }
  }

}


int main(int argc,  char *argv[])
{
  int	rank, nprocs, M=MAXSIZE;
  double tmpS[MAXSIZE];
  double tmpR[MAXSIZE];

  typedef Givaro::Modular<double> Field;
  Field F(17);

  MPI_Init(&argc,&argv);
  
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
 
  /*bool visited[nprocs];
  for(size_t j=0;j<nprocs;j++){
    visited[j] = false;
  }*/

  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  

  // Let A be a M times M random square matrix
  Field::Element_ptr A;
  A = fflas_new(F,M,1);
  RandomMatrix (F, M, 1, A, 1);

  MPI_Request request[nprocs];
  MPI_Status status[nprocs];
  int	prev;
  int	next;

//Determine the destination index for each rank
  prev = (rank == 0)?(nprocs-1):(rank-1);
  next = (rank == (nprocs - 1))?0:(rank+1);
//data generation by each process
    for(int count=0; count<MAXSIZE; count++){
      //tmpS[count] = rank;
      tmpR[count] = -1;
    }//buf[MAXSIZE]=rank;
unpack(F, A,  M, 1, 1, tmpS);
	std::cout<<"Process("<<rank<<") should send data to process("<<next<<"): "<<std::endl;
		    for(size_t j=0;j<MAXSIZE;j++){
		      std::cout<<""<<tmpS[j]<<"\n";
		    }std::cout<<std::endl;sleep(2); 


for(int step=0; step<nprocs-1;step++){ //loop controls the ring direction

//Nonblocking communication will send and receive data in disorder 
//MPI_Isend(&tmpS, MAXSIZE, MPI_DOUBLE, next, 123, MPI_COMM_WORLD, &request[rank]);  //MPI_Barrier(MPI_COMM_WORLD);
//MPI_Irecv(&tmpR, MAXSIZE, MPI_DOUBLE, prev, 123, MPI_COMM_WORLD, &request[prev]);//MPI_Wait(&request, &status);

	//MPI_Wait( &request[prev], NULL );

//Blocking communication will keed the send and receive order
MPI_Send(&tmpS, MAXSIZE, MPI_DOUBLE, next, 123, MPI_COMM_WORLD);
MPI_Recv(&tmpR, MAXSIZE, MPI_DOUBLE, prev, 123, MPI_COMM_WORLD, &status[rank]);//MPI_Wait(&request, &status);


	    //for(int i=0; i<rank-1; i++){
		 // for(int j=0; j<rank-1; j++){
		      //compute();
		   // }
	    //}	
repack(F, A, M, 1, 1, tmpR);
	std::cout<<"Process("<<rank<<") received data : step("<<step<<")"<<std::endl;
 
		    //for(size_t j=0;j<MAXSIZE;j++){
		      //std::cout<<""<<static_cast<int>(tmpR[j])<<"\t";
		  //  }std::cout<<std::endl;sleep(2);
FFLAS::WriteMatrix (std::cout << ""<< std::endl, F, M, 1, A, 1) << std::endl;


	double tmp;  
	//Swap(tmpS,tmpR);
	    for(int j=0; j<M; j++){
	      tmp=tmpS[j];
	      tmpS[j]=tmpR[j];
	      tmpR[j]=tmp;
	    }  
        sleep(2);
}

 
  
  MPI_Finalize();
  return 0;
}



