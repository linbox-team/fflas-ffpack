// ==========================================================================
// $Source: /var/lib/cvs/DLAFF/EXPERIMENTS/src/BLOCKING/mulMM.C,v $
// Copyright(c)'94-97 by Givaro Team
// Author: T. Gautier
// $Id: mulMM.C,v 1.1 2007-06-19 14:22:16 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
clock_t totaltime;


typedef long Base;
long P;
double dP;

inline void reduce(long& r) { r %= P; }
//inline void reduce(double& r) { r = fmod(r,dP); }


template<class T>
class Matrix {
    const int sq;
    const int size;
public:
  Matrix(int s): size(s),sq(s*s), _data(new T[sq]) {};
  int rowdim() const { return size; }
  int coldim() const { return size; }
  inline T& operator[](int i) { return _data[i]; }
  inline T& operator()(int i, int j) { return _data[i*size+j]; }
  inline const T& operator()(int i, int j) const { return _data[i*size+j]; }
  inline Matrix<T>& operator= (const T& val ) 
  { for (int i=size*size; --i; ) _data[i] =val; return *this; }
  inline void mul( const Matrix<T>& A, const Matrix<T>& B)
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

  inline void mulmod( const Matrix<T>& A, const Matrix<T>& B)
  {
      mul(A,B);
      for(int i=0; i<this->sq; ++i)
          reduce(this->_data[i]);
  }


  inline void addmul( const Matrix<T>& A, const Matrix<T>& B)
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
  inline void addmulmod( const Matrix<T>& A, const Matrix<T>& B)
  {
      addmul(A,B);
      for(int i=0; i<this->sq; ++i)
          reduce(this->_data[i]);
  }

private:
  T * _data; 
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
    for (i=0; i<DIM; i++) {
      for (j=0 ; j<DIM ; j++) {
        AA(i,j) = rand();
        BB(i,j) = rand();
        CC(i,j) = rand();
      }
    }

    CC.mul(AA,BB);
    totaltime=clock();
    for (i=0; i<NB; i++) CC.mul(AA,BB);
    totaltime=clock()-totaltime;
    seconds = (double)totaltime/CLOCKS_PER_SEC;

    std::cerr << "\nStandard Dim:" << DIM << ", Nb: " << NB << std::endl;
    std::cerr << "Mops:" << coef/seconds << std::endl;
    std::cerr << "time:" << seconds/NB << std::endl;
    MoyenneDesTemps += seconds;



    CC.mulmod(AA,BB);
    totaltime=clock();
    for (i=0; i<NB; i++) CC.mulmod(AA,BB);
    totaltime=clock()-totaltime;
    seconds = (double)totaltime/CLOCKS_PER_SEC;

    std::cerr << "\nStandardMod Dim:" << DIM << ", Nb: " << NB << std::endl;
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
