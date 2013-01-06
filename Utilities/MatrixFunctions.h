// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef UTILITIES_MATRIXFUNCTIONS_
#define UTILITIES_MATRIXFUNCTIONS_

#include<vector>
#include<cmath>
#include<stdio.h>
#include<string>
#include<iostream>
#include<fstream>
#include "Matrix.h"

// Utility functions and functions commonly used for matrices. Most functions
// overloaded to accept either a matrix object or a vector of vectors. Note that
// many of the functions that operate on matrix objects are slower than
// necessary because they simply extract a vector of vectors from the matrix and
// call the corresponding version of the code.

namespace utilities
{

// Prints out the matrix in a row/col format. The format string corresponds to
// the format in the printf function.
template <typename T> 
void PrintMatrix(const std::vector<std::vector<T> > &matrix, std::string format)
{
  for(unsigned int i = 0; i < matrix.size(); ++i)
  {
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
    {
      printf(format.c_str(), matrix[i][j]);
      printf(" ");
    }
    printf("\n");
  }
}

template <typename T> 
void PrintMatrix(const Matrix<T> &matrix, std::string format)
{
  std::vector<std::vector<T> > m = matrix.GetVectorOfVectors();
  PrintMatrix(m, format);
}

// Allows only a subset of the full matrix to be easily printed.
template <typename T> 
void PrintSubMatrix(const std::vector<std::vector<T> > &matrix, 
    unsigned int startr, unsigned startc, unsigned int rows, 
    unsigned int cols, std::string format)
{
  for(unsigned int i = startr; i < rows; ++i)
  {
    for(unsigned int j = startc; j < cols; ++j)
    {
      printf(format.c_str(), matrix[i][j]);
      printf(" ");
    }
    printf("\n");
  }
}

template <typename T> 
void PrintSubMatrix(const Matrix<T> &matrix, 
    unsigned int startr, unsigned startc, unsigned int rows, 
    unsigned int cols, std::string format)
{
  std::vector<std::vector<T> > m;
  m = matrix.GetVectorOfVectors();
  PrintSubMatrix(m, startr, startc, rows, cols, format);
}

template <typename T>
void MatrixAbs(std::vector<std::vector<T> > &matrix)
{
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      if( matrix[i][j] < 0)
        matrix[i][j] = -matrix[i][j];
}

template <typename T>
void MatrixAbs(Matrix<T> &matrix)
{
  for(unsigned int r = 0; r < matrix.NumRows(); ++r)
    for(unsigned int c = 0; c < matrix.NumCols(); ++c)
      matrix(r,c) = std::abs(matrix(r,c));
}

// Frobenius norm is the square root of the sums of the squares of each element.
// sqrt( sum( x^2) ).
template <typename T>
double FrobeniusNorm (const std::vector<std::vector<T> > &matrix)
{
  double ret = 0;
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      ret += static_cast<double>(matrix[i][j] * matrix[i][j]);

  return std::sqrt(ret);
}

template <typename T>
double FrobeniusNorm (const Matrix<T> &matrix)
{
  std::vector<std::vector<T> > m = matrix.GetVectorOfVectors();
  return FrobeniusNorm(m);
}

template <typename T>
void InitializeMatrix(std::vector<std::vector<T> > &matrix, T value, 
    unsigned int rows, unsigned int cols)
{
  matrix.clear();
  matrix.resize(rows);
  for(unsigned int i = 0; i < matrix.size(); ++i)
  {
    matrix[i].resize(cols, value);
  }
}

// Every element in the matrix is 0, except for the diagonal, which are 1.
template <typename T>
void CreateIdentityMatrix(std::vector<std::vector<T> > &matrix, 
    unsigned int rows, unsigned int cols)
{
  InitializeMatrix(matrix, static_cast<T>(0), rows, cols);
  unsigned int minval = rows;
  if( cols < rows)
    minval = cols;
  for(unsigned int i = 0; i < matrix.size(); ++i)
    matrix[i][i] = 1;
}

template<typename T>
void CreateIdentityMatrix(Matrix<T> &matrix, unsigned int rows, 
    unsigned int cols)
{
  matrix.Initialize(rows, cols, 0);
  matrix.SetDiagonal(1);
}

// Loops through the entire matrix to find the maximal element. Also gives the 
// location of the element.
template <typename T>
T MaxElementInMatrix(std::vector<std::vector<T> > &matrix, unsigned int &row, 
    unsigned int &col)
{
  T max_element = matrix[0][0];
  row = 0;
  col = 0;
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      if(matrix[i][j] > max_element)
      {
        max_element = matrix[i][j];
        row = i;
        col = j;
      }
  return max_element;
}

template <typename T>
T MaxElementInMatrix(const Matrix<T> &matrix, unsigned int &row, 
    unsigned int &col)
{
  std::vector<std::vector<T> > m = matrix.GetVectorOfVectors();
  return MaxElementInMatrix(m, row, col);
}

// Loops through the entire matrix to find the minimal element. Also gives the 
// location of the element.
template <typename T>
T MinElementInMatrix(std::vector<std::vector<T> > &matrix, unsigned int &row, 
    unsigned int &col)
{
  T min_element = matrix[0][0];
  row = 0;
  col = 0;
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      if(matrix[i][j] < min_element)
      {
        min_element = matrix[i][j];
        row = i;
        col = j;
      }
  return min_element;
}

template <typename T>
T MinElementInMatrix(const Matrix<T> &matrix, unsigned int &row, 
    unsigned int &col)
{
  std::vector<std::vector<T> > m = matrix.GetVectorOfVectors();
  return MinElementInMatrix(m, row, col);
}

// Standard matrix transposition.
template <typename T>
std::vector<std::vector<T> > TransposeMatrix(
    std::vector<std::vector<T> > &matrix)
{
  std::vector<std::vector<T> > ret;
  InitializeMatrix(ret, matrix[0][0], matrix.size(), matrix[0].size());
  for(unsigned int i = 0; i < matrix.size() - 1; ++i)
    for(unsigned int j = i+1; j < matrix[i].size(); ++j)
    {
      ret[i][j] = matrix[j][i];
      ret[j][i] = matrix[i][j];
    }
  return ret;
}

template <typename T>
Matrix<T> TransposeMatrix(const Matrix<T> &matrix)
{
  std::vector<std::vector<T> > ret, m;
  Matrix<T> matrix_ret;
  m = matrix.GetVectorOfVectors();
  ret = TransposeMatrix(m);
  matrix_ret.Initialize(ret);
  return matrix_ret;
}

// Uses a simple, but inefficient matrix product implementation. If this is
// something you need to do often with large matrices, you should be using a
// more powerful matrix library. Function does not check that the sizes of the
// matrices are valid for the product.
template <typename T>
std::vector<std::vector<T> > MatrixProduct(
    std::vector<std::vector<T> > &A, std::vector<std::vector<T> > &B)
{
  std::vector<std::vector<T> > ret;
  InitializeMatrix(ret, static_cast<T>(0), A.size(), B[0].size());
  for(unsigned int i = 0; i < A.size(); ++i)
    for(unsigned int j = 0; j < B[i].size(); ++j)
      for(unsigned int m = 0; m < B.size(); ++m)
        ret[i][j] += ( A[i][m] * B[m][j] );
  return ret;
}

template <typename T>
Matrix<T> MatrixProduct(const Matrix<T> &A, const Matrix<T> &B)
{
  std::vector<std::vector<T> > a = A.GetVectorOfVectors();
  std::vector<std::vector<T> > b = B.GetVectorOfVectors();
  std::vector<std::vector<T> > prod;
  Matrix<T> ret;
  prod = MatrixProduct(a, b);
  ret.Initialize(prod);
  return ret;
}

template <typename T>
std::vector<T> MatrixDiagonal(std::vector<std::vector<T> > &matrix)
{
  std::vector<T> ret;
  unsigned int min_dimension = matrix.size();
  if(matrix[0].size() < min_dimension)
    min_dimension = matrix[0].size();
  for(unsigned int i = 0; i < min_dimension; ++i)
    ret.push_back(matrix[i][i]);

  return ret;
}

template <typename T>
T SumMatrixValue(std::vector<std::vector<T> > &matrix)
{
  T ret = static_cast<T>(0);
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      ret += matrix[i][j];
  return ret;
}

template <typename T>
std::vector<std::vector<T> > UpperTriangularMatrix(
    std::vector<std::vector<T> > &matrix)
{
  std::vector<std::vector<T> > ret;
  InitializeMatrix(ret, static_cast<T>(0), matrix.size(), matrix[0].size());
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = i; j < matrix[i].size(); ++j)
      ret[i][j] = matrix[i][j];

  return ret;     
}

// Saves the matrix to a file in the row/col format. Each value is separate by a
// space and each row is given a newline.
template <typename T>
void SaveMatrixToTextFile(std::vector<std::vector<T> > &matrix, 
    std::string filename)
{
  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::out);
  for(unsigned int i = 0; i < matrix.size(); ++i)
  {
    fout<<matrix[i][0];
    for(unsigned int j = 1; j < matrix[i].size(); ++j)
      fout<<" "<<matrix[i][j];
    fout<<"\n";
  }

  fout.close();
}

template <typename T>
void SaveMatrixToTextFile(const Matrix<T> &matrix, std::string filename)
{
  std::vector<std::vector<T> > m = matrix.GetVectorOfVectors();
  SaveMatrixToTextFile(m, filename);
}

} // end namespace utilities
#endif
