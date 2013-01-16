#include<vector>
#include<iostream>
#include "Matrix.h"

void PrintMatrix(utilities::Matrix<char> m)
{
  for(unsigned int r = 0; r < m.NumRows(); ++r)
  {
    for(unsigned int c = 0; c < m.NumCols(); ++c)
      std::cout<<m(r,c)<<" ";
    std::cout<<"\n";
  }
}

int main()
{
  utilities::Matrix<char> a, b(3,4);
  a.Initialize(3, 2);
  b.Initialize(4, 4, 0);

  a(0,0) = 'a'; a(0,1) = 'b';
  a(1,0) = 'c'; a(1,1) = 'd';
  a(2,0) = 'e'; a(2,1) = 'f';

  PrintMatrix(a);
  a.Transpose();
  PrintMatrix(a);

  return 0;
}
