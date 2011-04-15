/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
#ifndef __MATIO_H
#define __MATIO_H
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
// Reading and writing matrices over double

// Reading a matrice from a (eventually zipped) file
double * read_dbl(char * mat_file,int* tni,int* tnj) {
  // char *UT;
  // int is_gzipped = 0;
  // size_t s = strlen(mat_file);
  double* X=NULL;
  const char * File_Name = mat_file;

  FILE* FileDes = fopen(File_Name, "r");
  if (FileDes != NULL) {
    char * tmp = new char[200];// usigned long tni, tnj;
    fscanf(FileDes,"%d %d %s\n",tni, tnj, tmp) ;
    int n=*tni;
    int p=*tnj;
    X = new double[n*p];
    for (int i=0;i<n*p;++i)
      X[i] = (double) 0;
    long i,j; long val;
    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    while(i && j) {
      X[p*(i-1)+j-1] = (double) val;
      fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    }
  }

  fclose(FileDes);
  return X;
}

// Displays a matrix
std::ostream& write_dbl(std::ostream& c,
		   double* E,
		   int n, int m, int id){

  for (int i = 0; i<n;++i){
    for (int j=0; j<m;++j)
      c << *(E+j+id*i) << " ";
    c << std::endl;
  }
  return c << std::endl;
}

// Reading and writing matrices over field

// Reading a matrice from a (eventually zipped) file
template<class Field>
typename Field::Element * read_field(const Field& F,char * mat_file,int* tni,int* tnj)
{
  // char *UT;
  // int is_gzipped = 0;
  // size_t s = strlen(mat_file);
  typename Field::Element zero;
  F.init(zero,0.0);
  typename Field::Element * X=NULL;
  const char * File_Name = mat_file;
  FILE* FileDes = fopen(File_Name, "r");
  if (FileDes != NULL) {
    char  tmp[200];// usigned long tni, tnj;
    fscanf(FileDes,"%d %d %s\n",tni, tnj, tmp) ;

    int n=*tni;
    int p=*tnj;
    X = new typename Field::Element[n*p];
    for (int i=0;i<n*p;++i)
      X[i] = zero;
    long i,j; long val;
    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    while(i && j) {
      F.init(X[p*(i-1)+j-1],double(val));
      fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    }
  }

  fclose(FileDes);
  return X;
}

// Displays a matrix
template<class Field>
std::ostream& write_field(const Field& F,std::ostream& c,
		     const typename Field::Element* E,
		     int n, int m, int id){

  typename Field::Element tmp;
  //#if DEBUG
  for (int i = 0; i<n;++i){
    for (int j=0; j<m;++j){
      F.convert(tmp,*(E+j+id*i));
      c << tmp << " ";}
    c << std::endl;
  }
  return c << std::endl;
  //#endif
}
#endif
