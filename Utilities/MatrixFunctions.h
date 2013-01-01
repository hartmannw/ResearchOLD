#ifndef UTILITIES_MATRIXFUNCTIONS_
#define UTILITIES_MATRIXFUNCTIONS_

#include<vector>
#include<cstdlib>
#include<cmath>
#include<stdio.h>
#include<string>
#include<iostream>
#include<fstream>

namespace utilities
{

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
void MatrixAbs(std::vector<std::vector<T> > &matrix)
{
  for(unsigned int i = 0; i < matrix.size(); ++i)
    for(unsigned int j = 0; j < matrix[i].size(); ++j)
      if( matrix[i][j] < 0)
        matrix[i][j] = -matrix[i][j];
}

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
}
#endif
