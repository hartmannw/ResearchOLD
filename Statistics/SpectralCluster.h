#ifndef STATISTICS_SPECTRALCLUSTER_
#define STATISTICS_SPECTRALCLUSTER_

#include<vector>
#include<cmath>
#include<algorithm>
#include<iostream>
#include "MatrixFunctions.h"

namespace statistics
{

class SpectralCluster
{
 private:
  std::vector<std::vector<double> > similarity_matrix_;
  static const double EPS_ = 10e-08;

  bool Jacobi(std::vector<std::vector<double> > &vectors, 
      std::vector<double> &values);
  bool Jacobi2(std::vector<std::vector<double> > &U, std::vector<double> &e);
  void ZeroRowCol(std::vector<std::vector<double> > &S,
      std::vector<std::vector<double> > &vectors, unsigned int row, 
      unsigned int col);
  bool Schur(std::vector<std::vector<double> > &U, std::vector<double> &e);
  bool FindNextPivot(std::vector<std::vector<double> > &M, 
      unsigned int &row, unsigned int &col);
 
 public:
  SpectralCluster(){}
  ~SpectralCluster(){}

  void SetMatrix(std::vector<std::vector<double> > matrix){
      similarity_matrix_ = matrix;}
  void ConvertToSimpleEigenvalue();

  std::vector<std::vector<double> > similarity_matrix(){
      return similarity_matrix_;}

  bool PerformEigenAnalysis(std::vector<std::vector<double> > &vectors, 
      std::vector<double> &values);

  void SortEigenvalues(std::vector<std::vector<double> > &vectors,
      std::vector<double> &values);

};

}

#endif
