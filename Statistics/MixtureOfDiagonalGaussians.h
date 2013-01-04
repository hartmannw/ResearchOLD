#ifndef STATISTICS_MIXTUREOFDIAGONALGAUSSIANS_
#define STATISTICS_MIXTUREOFDIAGONALGAUSSIANS_

#include <vector>
#include <cmath>
#include <iostream>
#include "DiagonalGaussian.h"

namespace statistics
{

class MixtureOfDiagonalGaussians
{
 private:
  std::vector<DiagonalGaussian> gaussian_;
  std::vector<double> weight_;
  
 public:
  MixtureOfDiagonalGaussians(){}
  ~MixtureOfDiagonalGaussians(){}

  void AddGaussian(DiagonalGaussian g, double w){
      gaussian_.push_back(g); weight_.push_back(w);}
  void SetAllGaussians(std::vector<DiagonalGaussian> gaussian, 
      std::vector<double> weight);
  double Likelihood(std::vector<double> &point);
  double LogLikelihood(std::vector<double> &point);

  double CSDivergence(MixtureOfDiagonalGaussians mog);

  unsigned int components(){return gaussian_.size();}
  double weight(unsigned int i){return weight_[i];}
  DiagonalGaussian gaussian(unsigned int i){return gaussian_[i];}
 
};

}

#endif
