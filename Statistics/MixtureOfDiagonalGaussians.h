// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_MIXTUREOFDIAGONALGAUSSIANS_
#define STATISTICS_MIXTUREOFDIAGONALGAUSSIANS_

#include <vector>
#include <cmath>
#include <iostream>
#include "DiagonalGaussian.h"

// Stores and evaluates mixtures of Gaussians where the covariance matrix is
// diagonal. This code does not have the ability to learn the parameters of the
// Gaussians; we assume the Gaussians have already been trained.

namespace statistics
{

class MixtureOfDiagonalGaussians
{
 private:
  std::vector<DiagonalGaussian> gaussian_;
  std::vector<double> weight_; // Weight vector should be normalized to one, but
                               // this is not enforced.
  
 public:
  MixtureOfDiagonalGaussians(){}
  ~MixtureOfDiagonalGaussians(){}

  // Adds an additional Gaussian to the mixture. The weight must be included,
  // but no constraints are enforced on the value.
  void AddGaussian(DiagonalGaussian g, double w){
      gaussian_.push_back(g); weight_.push_back(w);}

  // Initializes the mixture by supplying a vector of Gaussians and the weight
  // vector.
  void SetAllGaussians(std::vector<DiagonalGaussian> gaussian, 
      std::vector<double> weight);

  // Returns the likelihood of the entire mixture.
  double Likelihood(std::vector<double> &point);
  double LogLikelihood(std::vector<double> &point);

  // Returns a value known as the Cauchy-Schwarz divergence, a symmetric
  // distance measure between two MOG distributions. It is similar to KL 
  // divergence, but it has a closed form solution and can be computed quickly.
  // More information can be found in the paper "Closed-form Cauchy-Schwarz PDF
  // divergence for mixture of Gaussians" by K. Kampa et al.
  double CSDivergence(MixtureOfDiagonalGaussians mog);

  // Returns number of Gaussians in the mixture.
  unsigned int components(){return gaussian_.size();}

  // Standard accessor functions.
  double weight(unsigned int i){return weight_[i];}
  DiagonalGaussian gaussian(unsigned int i){return gaussian_[i];}
 
};

}

#endif
