// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef UTILITIES_MATRIX_H_
#define UTILITIES_MATRIX_H_

#include<iostream>
#include<vector>
#include<algorithm>

namespace utilities
{

// Matrix is used for storing data in a two-dimensional vector of a size defined
// at construction. It is designed to replace the vector of vectors approach 
// when the matrix is a static size. 
template <class T>
class Matrix
{
 private:
  std::vector<T> matrix_;
  unsigned int rows_, cols_;

 public:
  Matrix() : rows_(0), cols_(0) {}
  Matrix(unsigned int rows, unsigned int cols) : matrix_(rows*cols), 
      rows_(rows), cols_(cols)  {}
  Matrix(unsigned int rows, unsigned int cols, T value);
  Matrix(const std::vector<std::vector<T> > &matrix);
  ~Matrix() {}

  bool Initialize(unsigned int rows, unsigned int cols);
  bool Initialize(unsigned int rows, unsigned int cols, T value);
  bool Initialize(const std::vector<std::vector<T> > &matrix);

  T operator() (unsigned int row, unsigned int col) const;
  T& operator() (unsigned int row, unsigned int col);
  unsigned int NumRows() const { return rows_;}
  unsigned int NumCols() const { return cols_;}
  std::vector<T> GetRow(unsigned int row) const;
  std::vector<T> GetCol(unsigned int col) const;
  std::vector<T> GetDiagonal() const;
  std::vector<std::vector<T> > GetVectorOfVectors() const;

  bool SetRow(unsigned int row, const std::vector<T> &values);
  bool SetRow(unsigned int row, T value);
  bool SetCol(unsigned int col, const std::vector<T> &values);
  bool SetCol(unsigned int col, T value);
  bool SetDiagonal(const std::vector<T> &values);
  bool SetDiagonal(T value);

  bool isSquare() const { return rows_ == cols_; }
  bool Transpose();
};

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, T value)
{
  Initialize(rows, cols, value);
}

template<class T>
Matrix<T>::Matrix(const std::vector<std::vector<T> > &matrix)
{
  Initialize(matrix);
}

template<class T>
bool Matrix<T>::Initialize(unsigned int rows, unsigned int cols)
{
  rows_ = rows;
  cols_ = cols;
  matrix_.resize(rows * cols);
  return true;
}

template<class T>
bool Matrix<T>::Initialize(unsigned int rows, unsigned int cols, T value)
{
  rows_ = rows;
  cols_ = cols;
  matrix_.clear();

  matrix_.resize(rows * cols);
  for(unsigned int i = 0; i < matrix_.size(); ++i)
    matrix_[i] = value;
  return true;
}

template<class T>
bool Matrix<T>::Initialize(const std::vector<std::vector<T> > &matrix)
{
  rows_ = matrix.size();
  if( rows_ == 0)
    return Initialize(0, 0);
  cols_ = matrix[0].size();
  matrix_.resize(rows_ * cols_);
  std::cout<<rows_<<" "<<cols_<<std::endl;
  for(unsigned int r = 0; r < rows_; ++r)
  {
    if(matrix[r].size() != cols_) // Each row of the matrix does not have the
    {                             // same number of columns.
      Initialize(0, 0);
      return false;
    }
    for(unsigned int c = 0; c < cols_; ++c)
    {
      matrix_[ (cols_ * r) + c ] = matrix[r][c];
    }
  }
  return true;
}

template<class T>
T Matrix<T>::operator() (unsigned int row, unsigned int col) const
{
  return matrix_[ (row * cols_) + col ];
}

template<class T>
T& Matrix<T>::operator() (unsigned int row, unsigned int col)
{
  return matrix_[ (row * cols_) + col ];
}

template<class T>
std::vector<T> Matrix<T>::GetRow(unsigned int row) const
{
  std::vector<T> ret;
  if(row >= rows_) // This row does not exist;
    return ret;
  ret.resize(cols_);
  for(unsigned int c = 0; c < cols_; ++c)
    ret[c] = matrix_[ (row * cols_) + c ];
  return ret;
}

template<class T>
std::vector<T> Matrix<T>::GetCol(unsigned int col) const
{
  std::vector<T> ret;
  if(col >= cols_) // This col does not exist;
    return ret;
  ret.resize(rows_);
  for(unsigned int r = 0; r < rows_; ++r)
    ret[r] = matrix_[ (r * cols_) + col ];
  return ret;
}

template<class T>
std::vector<std::vector<T> > Matrix<T>::GetVectorOfVectors() const
{
  std::vector<std::vector<T> > ret;
  ret.resize(rows_);
  for(unsigned int r = 0; r < rows_; ++r)
  {
    ret[r].resize(cols_);
    for(unsigned int c = 0; c < cols_; ++c)
      ret[r][c] = matrix_[ (r * cols_) + c];
  }
  return ret;
}

template<class T>
bool Matrix<T>::SetRow(unsigned int row, const std::vector<T> &values)
{
  if(row > rows_) // This row does not exist
    return false;
  if(values.size() != cols_) // Vector does not match number of cols
    return false;
  for(unsigned int c = 0; c < cols_; ++c)
    matrix_[ (row * cols_) + c ] = values[c];
  return true;
}

template<class T>
bool Matrix<T>::SetRow(unsigned int row, T value)
{
  if(row > rows_) // This row does not exist
    return false;
  for(unsigned int c = 0; c < cols_; ++c)
    matrix_[ (row * cols_) + c ] = value;
  return true;
}

template<class T>
bool Matrix<T>::SetCol(unsigned int col, const std::vector<T> &values)
{
  if(col > cols_) // This col does not exist
    return false;
  if(values.size() != rows_) // Vector does not match number of rows
    return false;
  for(unsigned int r = 0; r < rows_; ++r)
    matrix_[ (r * cols_) + col ] = values[r];
  return true;
}

template<class T>
bool Matrix<T>::SetCol(unsigned int col, T value)
{
  if(col > cols_) // This row does not exist
    return false;
  for(unsigned int r = 0; r < rows_; ++r)
    matrix_[ (r * cols_) + col ] = value;
  return true;
}

template<class T>
bool Matrix<T>::SetDiagonal(const std::vector<T> &values)
{
  if( ! isSquare() ) // Cannot set the diagonal of a nonsquare matrix.
    return false;
  if(values.size() != rows_) // Vector length does not match diagonal size.
    return false;
  for(unsigned int i = 0; i < rows_; ++i)
    matrix_[ (i * cols_) + i] = values[i];
  return true;
}

template<class T>
bool Matrix<T>::SetDiagonal(T value)
{
  if( ! isSquare() ) // Cannot set the diagonal of a nonsquare matrix.
    return false;
  for(unsigned int i = 0; i < rows_; ++i)
    matrix_[ (i * cols_) + i] = value;
  return true;
}

// Note that a transposition is really just a representation issue. It is 
// possible and probably more efficient to simply set a transposition flag. Then
// every interaction with the class flips the rows and cols. It does add extra
// to every function call, so for the time being we are performing the less 
// efficient action of actually moving the data.
template<class T>
bool Matrix<T>::Transpose()
{
  std::vector<T> new_matrix = matrix_;
  for(unsigned int r = 0; r < rows_; ++r)
    for(unsigned int c = 0; c < cols_; ++c)
    {
      unsigned int index = (r * cols_) + c;
      unsigned int transpose = (c * rows_) + r;
      new_matrix[transpose] = matrix_[index];
    }
  std::swap(rows_, cols_);
  matrix_ = new_matrix;
  return true;
}

} // end namespace utilities
#endif
