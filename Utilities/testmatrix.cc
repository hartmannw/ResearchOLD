#include<vector>
#include<iostream>
#include "Matrix.h"

int main()
{
  utilities::Matrix<double> a, b(3,4);
  a.Initialize(3, 3);
  b.Initialize(4, 4, 0);

  std::cout<<b(2,2)<<std::endl;
  b(2,2) = 3.1415926;
  std::cout<<b(2,2)<<std::endl;
  return 0;
}
