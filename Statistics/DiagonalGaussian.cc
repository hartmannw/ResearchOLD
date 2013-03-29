// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "DiagonalGaussian.h"

namespace statistics
{

void DiagonalGaussian::Initialize(std::vector<double> mean, 
    std::vector<double> variance)
{
  mean_ = mean;
  variance_ = variance;
  CalculateConstant();
}

void DiagonalGaussian::CalculateConstant()
{
  double pi = 4.0 * atan(1.0); // Standard method of calculating pi.
  constant_ = 1;
  for(unsigned int i =0; i < variance_.size(); ++i)
    constant_ *= variance_[i];
  constant_ = 1 / (std::sqrt(constant_) * std::pow(2*pi, mean_.size() / 2.0));
}

void DiagonalGaussian::AddMean(std::vector<double> mean)
{
  for(unsigned int i = 0; i < mean.size(); ++i)
    mean_[i] += mean[i];
  CalculateConstant();
}

void DiagonalGaussian::AddVariance(std::vector<double> variance)
{
  for(unsigned int i = 0; i < variance.size(); ++i)
    variance_[i] += variance[i];
  CalculateConstant();
}

double DiagonalGaussian::Likelihood(std::vector<double> point) const
{
  double result = 0;
  for(unsigned int i =0; i < point.size(); ++i)
    result += pow(mean_[i] - point[i],2.0) / (variance_[i]);

  result = constant_ *  std::exp(-0.5 * result);
  return result;
}

double DiagonalGaussian::LogLikelihood(std::vector<double> point) const
{
  return std::log(Likelihood(point));
}

double DiagonalGaussian::KLDivergence(DiagonalGaussian g)
{
  double ret = 0;
  double det_term = 1;
  for(unsigned int i = 0; i < mean_.size(); ++i)
  {
    ret += variance_[i] / g.variance(i);
    ret += pow(g.mean(i) - mean_[i], 2.0) / g.variance(i);
    det_term *= (variance_[i] / g.variance(i));
  }
  ret = 0.5 * (ret - std::log(det_term) - dimension());
  return ret;
}

double DiagonalGaussian::SymmetricKLDivergence(DiagonalGaussian g)
{
  return KLDivergence(g) + g.KLDivergence(*this);
}

// Determinant of a diagonal matrix is the product of the diagaonal elements.
double DiagonalGaussian::determinant() const
{
  double ret = 1;
  for(unsigned int i = 0; i < variance_.size(); ++i)
    ret *= variance_[i];
  return ret;
}


}
