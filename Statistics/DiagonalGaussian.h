// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_DIAGONALGUASSIAN_
#define STATISTICS_DIAGONALGUASSIAN_

#include <vector>
#include <cmath>

// Stores the information needed for a Gaussian with diagonal variance. Both the
// likelihood and log likelihood of the Gaussian can be evaluated. KL divergence
// and symmetric KL divergence are also implemented.
namespace statistics
{

class DiagonalGaussian
{
 public:
  DiagonalGaussian(){}
  ~DiagonalGaussian(){}

  // Creates the Gaussian given a mean and variance vector.
  void Initialize(std::vector<double> mean, std::vector<double> variance);

  // Calculates the normalization constant used determining the likelihood.
  void CalculateConstant();

  // Adds a vector of values to the mean vector of the Gaussian.
  void AddMean(std::vector<double> mean);

  // Adds a vector of values to the variance vector of the Gaussian.
  void AddVariance(std::vector<double> variance);

  // Returns the likelihood of the point given the Gaussian. The dimension of
  // point must be equally to the dimension of the Gaussian.
  double Likelihood(std::vector<double> point);
  double LogLikelihood(std::vector<double> point);

  // KL Divergence is a measure of the distance between two Gaussian
  // distributions.
  double KLDivergence(DiagonalGaussian g);          // I believe these functions
  double SymmetricKLDivergence(DiagonalGaussian g); // are correct, but they
                                                    // should be tested more.

  // Standard accessor functions.
  unsigned int dimension(){ return mean_.size();}
  double mean(unsigned int i){ return mean_[i];}
  double variance(unsigned int i){return variance_[i];}
  std::vector<double> mean(){return mean_;}
  std::vector<double> variance(){return variance_;}
  
  // Returns the determinant of the diagonal covariance matrix.
  double determinant();

 private:
  std::vector<double> mean_;
  std::vector<double> variance_; // Since the Gaussian has a digaonal covariance
                                 // matrix, it can be stored as a vector.
  double constant_; // standard Gaussian normalization term
};

}

#endif
