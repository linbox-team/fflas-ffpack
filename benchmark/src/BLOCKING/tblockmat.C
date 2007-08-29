// ==========================================================================
// $Source: /var/lib/cvs/DLAFF/EXPERIMENTS/src/BLOCKING/tblockmat.C,v $
// Copyright(c)'94-97 by Givaro Team
// Author: T. Gautier
// $Id: tblockmat.C,v 1.1 2007-06-19 14:22:16 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
clock_t totaltime;


// Taille des blocs
#ifndef NDIM
#define NDIM TOREPLACE
#endif


typedef long Base;
long P;
double dP;

inline void reduce(long& r) { r %= P; }
//inline void reduce(double& r) { r = fmod(r,dP); }


template<class T, int size>
class FixedMatrix {
    const int sq;
public:
  FixedMatrix():sq(size*size) {};
  int rowdim() const { return size; }
  int coldim() const { return size; }
  inline T& operator[](int i) { return _data[i]; }
  inline T& operator()(int i, int j) { return _data[i*size+j]; }
  inline const T& operator()(int i, int j) const { return _data[i*size+j]; }
  inline FixedMatrix<T,size>& operator= (const T& val ) 
  { for (int i=size*size; --i; ) _data[i] =val; return *this; }
  inline void mul( const FixedMatrix<T,size>& A, const FixedMatrix<T,size>& B)
  {
    register const T* Ai = A._data;
    T* Ci = _data;
    for (int i=size; --i; Ai += size, Ci += size)
    {
      for (int j=size; --j; )
      {
        register const T* Bj = &B._data[j];
        T sum =0;
        for (int k=size; --k; Bj+=size)
          sum += Ai[k] * *Bj;
        Ci[j] = sum;
      }
    }
  }
  inline void addmul( const FixedMatrix<T,size>& A, const FixedMatrix<T,size>& B)
  {
    register const T* Ai = A._data;
    T* Ci = _data;
    for (int i=size; --i; Ai += size, Ci += size)
    {
      for (int j=size; --j; )
      {
        register const T* Bj = &B._data[j];
        T sum =0;
        for (int k=size; --k; Bj+=size)
          sum += Ai[k] * *Bj;
        Ci[j] += sum;
      }
    }
  }
  inline void addmulmod( const FixedMatrix<T,size>& A, const FixedMatrix<T,size>& B)
  {
      addmul(A,B);
      for(int i=0; i<this->sq; ++i)
          reduce(this->_data[i]);
  }

private:
  T _data[size*size]; 
};

    

template<class T>    
class Matrix {                                                                            
public:                                                                                        
  Matrix(int dim) {
   _dim = dim;
   _nblock = dim / NDIM;
   _data = new FixedMatrix<T,NDIM>[_nblock*_nblock];
  } 
  ~Matrix() { delete [] _data; }
  int rowdim() const { return _dim; }
  int coldim() const { return _dim; }
  inline FixedMatrix<T,NDIM>& operator[](int i)
  { return _data[i]; }
  inline FixedMatrix<T,NDIM>& operator()(int i, int j) 
  { return _data[i*_nblock+j]; }
  inline const FixedMatrix<T,NDIM>& operator()(int i, int j) const 
  { return _data[i*_nblock+j]; }
  inline void mul( const Matrix<T>& A, const Matrix<T>& B)                 
  {
    for (int i=0; i<_nblock; ++i )                                              
      for (int j=0; j<_nblock; ++j )
      {
        (*this)(i,j)=0;
        for (int k=0; k<_nblock; ++k )  {
          (*this)(i,j).addmul( A(i,k), B(k,j));
        }
      }
  }
  inline void mulmod( const Matrix<T>& A, const Matrix<T>& B)                 
  {
    for (int i=0; i<_nblock; ++i )                                              
      for (int j=0; j<_nblock; ++j )
      {
        (*this)(i,j)=0;
        for (int k=0; k<_nblock; ++k )  {
          (*this)(i,j).addmulmod( A(i,k), B(k,j));
        }
      }
  }
private:
  int _dim;
  int _nblock;
  FixedMatrix<T,NDIM>* _data;
};                                                                                             


template<class T>
void print_mat( int dim, const T* A)
{
  int i,j;
  for (i=0; i<dim; i++)
  {
    for (j=0; j<dim; j++)
      std::cout << A[i+j*dim] << "  ";
    std::cout << std::endl;
  }
}


int main(int argc, char** argv)
{
  int DIM = atoi(argv[1]) ;
  int NB = (argc>2?atoi(argv[2]):1) ;
  P = (argc>3?atoi(argv[3]):65521);
  dP = (double)P;
  char * nomfich = new char[6+strlen(argv[0])]; sprintf(nomfich, "plot1-%s", argv[0]);
  std::cerr << nomfich << std::endl;
  std::ofstream plot1(nomfich, std::ios::app);
  delete [] nomfich;

  Matrix<Base> AA(DIM);
  Matrix<Base> BB(DIM);
  Matrix<Base> CC(DIM);
  int i,j ;
  srand((int)clock());
  int MAXITER = 3;
  double MoyenneDesTemps = 0.0;
  double MoyenneDesTempsMod = 0.0;
  double coef = DIM;
  coef = coef*coef*(2.0*coef-1)*1e-6*NB;
  double seconds;

  for (int k=0; k<MAXITER; k++) {  
    for (i=0; i<(DIM*DIM)/NDIM/NDIM; i++) {
      for (j=0 ; j<NDIM*NDIM ; j++) {
        AA[i][j] = rand();
        BB[i][j] = rand();
        CC[i][j] = rand();
      }
    }

    CC.mul(AA,BB);
    totaltime=clock();
    for (i=0; i<NB; i++) CC.mul(AA,BB);
    totaltime=clock()-totaltime;
    seconds = (double)totaltime/CLOCKS_PER_SEC;

    std::cerr << "\nMUL3 Dim:" << DIM << ", Nb: " << NB << std::endl;
    std::cerr << "Mops:" << coef/seconds << std::endl;
    std::cerr << "time:" << seconds/NB << std::endl;
    MoyenneDesTemps += seconds;



    CC.mulmod(AA,BB);
    totaltime=clock();
    for (i=0; i<NB; i++) CC.mulmod(AA,BB);
    totaltime=clock()-totaltime;
    seconds = (double)totaltime/CLOCKS_PER_SEC;

    std::cerr << "\nMULMOD3 Dim:" << DIM << ", Nb: " << NB << std::endl;
    std::cerr << "Mops:" << coef/seconds << std::endl;
    std::cerr << "time:" << seconds/NB << std::endl;
    MoyenneDesTempsMod  += seconds;
  }

  MoyenneDesTemps /= (double)(MAXITER);
  MoyenneDesTempsMod /= (double)(MAXITER);
  

  plot1 << DIM << '\t' << MoyenneDesTemps/(double)NB << '\t' << coef/MoyenneDesTemps;
  plot1 << '\t' << MoyenneDesTempsMod/(double)NB << '\t' << coef/MoyenneDesTempsMod << std::endl;
  plot1.close();
  return 0 ;
};
